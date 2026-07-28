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
#include "bu/time.h"
#include "bu/vls.h"
#include "raytrace.h"

#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <mutex>
#include <shared_mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#define POP_MAXLEVEL 16
#define POP_CACHEDIR BOBOL_DRAW_CACHE_DIR
#define CACHE_CURRENT_FORMAT 10

#define CACHE_POP_MAX_LEVEL "max"
#define CACHE_POP_MIN_LEVEL "min"
#define CACHE_VERTEX_COUNT "vc"
#define CACHE_TRI_COUNT "tc"
#define CACHE_OBJ_BOUNDS "bb"
#define CACHE_SHADED_CULL_BACKFACES "cb"
#define CACHE_VERT_LEVEL "v"
#define CACHE_VERTNORM_LEVEL "vn"
#define CACHE_TRI_LEVEL "t"

struct BObolMeshLodContextInternal {
    struct bu_cache *lodCache;
    struct bu_cache *nameCache;
    struct bu_vls *fname;
    char *registryKey;
    std::shared_mutex *accessMutex;
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
    const int *faces;
    int pcnt;
    const point_t *points;
    int porig_cnt;
    const point_t *points_orig;
    const vect_t *normals;
    point_t bmin;
    point_t bmax;
};

struct BObolMeshLodBotDetail {
    struct db_i *dbip;
    struct directory *dp;
    struct rt_db_internal *intern;
};

class BObolPopRec
{
public:
    unsigned short x = 0;
    unsigned short y = 0;
    unsigned short z = 0;
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

static std::map<std::string, struct BObolMeshLodContext *> &
mesh_lod_context_registry(void)
{
    static std::map<std::string, struct BObolMeshLodContext *> registry;
    return registry;
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

/* A consistently wound closed mesh may still contain independently oriented
 * shells.  For an explicitly unoriented BoT there is no material-side
 * contract telling us how those shells relate, so automatic exterior
 * normalization is safe only for one connected component.  This union/find
 * is linear and substantially smaller than the half-edge topology check that
 * precedes it. */
static bool
mesh_lod_faces_one_component(const int *faces, size_t faceCount,
			     size_t vertexCount)
{
    if (!faces || !faceCount || !vertexCount ||
	vertexCount > static_cast<size_t>(INT_MAX))
	return false;

    const int inactive = std::numeric_limits<int>::min();
    std::vector<int> parent(vertexCount, inactive);
    size_t components = 0;
    const auto activate = [&parent, &components, inactive](int vertex) {
	if (parent[static_cast<size_t>(vertex)] == inactive) {
	    parent[static_cast<size_t>(vertex)] = -1;
	    components++;
	}
    };
    const auto findRoot = [&parent](int vertex) {
	int root = vertex;
	while (parent[static_cast<size_t>(root)] >= 0)
	    root = parent[static_cast<size_t>(root)];
	while (vertex != root) {
	    const int next = parent[static_cast<size_t>(vertex)];
	    parent[static_cast<size_t>(vertex)] = root;
	    vertex = next;
	}
	return root;
    };
    const auto unite = [&parent, &components, &findRoot](int a, int b) {
	int rootA = findRoot(a);
	int rootB = findRoot(b);
	if (rootA == rootB)
	    return;
	/* Negative roots store component size. */
	if (parent[static_cast<size_t>(rootA)] >
	    parent[static_cast<size_t>(rootB)])
	    std::swap(rootA, rootB);
	parent[static_cast<size_t>(rootA)] +=
	    parent[static_cast<size_t>(rootB)];
	parent[static_cast<size_t>(rootB)] = rootA;
	components--;
    };

    for (size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex) {
	const int *face = &faces[faceIndex * 3];
	activate(face[0]);
	activate(face[1]);
	activate(face[2]);
	unite(face[0], face[1]);
	unite(face[0], face[2]);
    }
    return components == 1;
}

/* Return +1 for exterior CCW, -1 for inward winding, and zero when numerical
 * cancellation makes the orientation indeterminate.  Translating about one
 * mesh vertex keeps the determinant stable for models far from the origin. */
static int
mesh_lod_closed_winding_sign(const point_t *vertices, size_t vertexCount,
			     const int *faces, size_t faceCount)
{
    if (!vertices || !vertexCount || !faces || !faceCount)
	return 0;

    const point_t &origin = vertices[static_cast<size_t>(faces[0])];
    long double signedSixVolume = 0.0L;
    long double absoluteTerms = 0.0L;
    for (size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex) {
	const int *face = &faces[faceIndex * 3];
	long double a[3];
	long double b[3];
	long double c[3];
	for (int axis = 0; axis < 3; ++axis) {
	    a[axis] = static_cast<long double>(
		vertices[static_cast<size_t>(face[0])][axis] - origin[axis]);
	    b[axis] = static_cast<long double>(
		vertices[static_cast<size_t>(face[1])][axis] - origin[axis]);
	    c[axis] = static_cast<long double>(
		vertices[static_cast<size_t>(face[2])][axis] - origin[axis]);
	}
	const long double term =
	    a[0] * (b[1] * c[2] - b[2] * c[1]) -
	    a[1] * (b[0] * c[2] - b[2] * c[0]) +
	    a[2] * (b[0] * c[1] - b[1] * c[0]);
	signedSixVolume += term;
	absoluteTerms += std::fabs(term);
    }
    const long double tolerance =
	std::numeric_limits<long double>::epsilon() *
	std::max(1.0L, absoluteTerms) * 64.0L;
    if (std::fabs(signedSixVolume) <= tolerance)
	return 0;
    return signedSixVolume > 0.0L ? 1 : -1;
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
    internal->accessMutex = new (std::nothrow) std::shared_mutex;
    if (!internal->accessMutex) {
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

    internal->lodCache = bu_cache_open(bu_vls_cstr(&lodCachePath), 1, 0);
    internal->nameCache = bu_cache_open(bu_vls_cstr(&nameCachePath), 1, 0);

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

    return mesh_lod_context_create(ctxName);
}

static unsigned long long
mesh_lod_key_get(struct BObolMeshLodContext *context, const char *name)
{
    if (!context || !name)
	return 0;

    {
	std::shared_lock<std::shared_mutex> lock(
	    *context->i->accessMutex);
	const auto found = context->i->nameKeys->find(name);
	if (found != context->i->nameKeys->end())
	    return found->second;
    }

    struct bu_vls keystr = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&keystr, "%s", name);
    if (bu_vls_strlen(&keystr) < 10)
	bu_vls_printf(&keystr, "GGGGGGGGGGGGG");
    unsigned long long hash = bu_data_hash(bu_vls_cstr(&keystr),
					   bu_vls_strlen(&keystr) * sizeof(char));
    bu_vls_sprintf(&keystr, "%llu:namekey", hash);

    void *data = NULL;
    std::shared_lock<std::shared_mutex> lock(*context->i->accessMutex);
    size_t dsize = bu_cache_get(&data, bu_vls_cstr(&keystr),
				context->i->nameCache, NULL);
    lock.unlock();
    bu_vls_free(&keystr);
    if (!dsize || !data)
	return 0;

    unsigned long long meshKey = *(unsigned long long *)data;
    bu_free(data, "lod key data");
    {
	std::unique_lock<std::shared_mutex> memoryLock(
	    *context->i->accessMutex);
	(*context->i->nameKeys)[name] = meshKey;
    }
    return meshKey;
}

static int
mesh_lod_key_put(struct BObolMeshLodContext *context,
		 const char *name,
		 unsigned long long key)
{
    if (!context || !name || !key)
	return -1;

    struct bu_vls keystr = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&keystr, "%s", name);
    if (bu_vls_strlen(&keystr) < 10)
	bu_vls_printf(&keystr, "GGGGGGGGGGGGG");
    unsigned long long hash = bu_data_hash(bu_vls_cstr(&keystr),
					   bu_vls_strlen(&keystr) * sizeof(char));
    bu_vls_sprintf(&keystr, "%llu:namekey", hash);

    std::unique_lock<std::shared_mutex> lock(*context->i->accessMutex);
    size_t wsize = bu_cache_write((void *)&key, sizeof(key),
				  bu_vls_cstr(&keystr),
				  context->i->nameCache, NULL);
    if (wsize > 0)
	(*context->i->nameKeys)[name] = key;
    lock.unlock();
    bu_vls_free(&keystr);
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
    ~BObolPopState();

    int getLevel(fastf_t len);
    bool setLevel(int level, bool materializeSnapped = true);
    void shrinkMemory(void);
    void levelPoint(point_t *out, const point_t *point, int level);
    void hierarchyInfo(struct BObolMeshLodHierarchyInfo *info) const;

    std::vector<int> lodTris;
    std::vector<fastf_t> lodTriPoints;
    std::vector<fastf_t> lodTriPointsSnapped;
    std::vector<fastf_t> lodTriNormals;

    int currLevel = -1;
    bool forceUpdate = false;
    int minPopLevel = 0;
    int maxPopLevel = 0;
    bool isValid = false;
    bool shadedCullBackfaces = false;
    unsigned long long hash = 0;
    point_t bbmin = VINIT_ZERO;
    point_t bbmax = VINIT_ZERO;
    struct BObolMeshLod *lod = NULL;

    BObolMeshLodDetailSetupCallback fullDetailSetup = NULL;
    BObolMeshLodDetailClearCallback fullDetailClear = NULL;
    BObolMeshLodDetailFreeCallback fullDetailFree = NULL;
    void *detailData = NULL;

private:
    void triProcess(void);
    int toLevel(int value, int level);
    fastf_t snap(fastf_t value, fastf_t min, fastf_t max, int level);
    bool isEqual(BObolPopRec r1, BObolPopRec r2, int level);
    bool triDegenerate(BObolPopRec r0, BObolPopRec r1,
		       BObolPopRec r2, int level);
    void cache(void);
    bool cacheTri(void);
    bool cacheWrite(const char *component, std::stringstream &stream);
    bool cacheWriteData(const char *component, const void *data, size_t size);
    size_t cacheGet(void **data, const char *component);
    void cacheDone(void);
    bool triPopLoad(int startLevel, int level, bool materializeSnapped);
    void triPopTrim(int level, bool materializeSnapped);
    void updateSnappedPoints(int level);
    bool setupFullDetail(void);
    void clearFullDetail(void);

    fastf_t minx = std::numeric_limits<fastf_t>::max();
    fastf_t miny = std::numeric_limits<fastf_t>::max();
    fastf_t minz = std::numeric_limits<fastf_t>::max();
    fastf_t maxx = -std::numeric_limits<fastf_t>::max();
    fastf_t maxy = -std::numeric_limits<fastf_t>::max();
    fastf_t maxz = -std::numeric_limits<fastf_t>::max();
    std::vector<unsigned short> precomputedMasks;
    size_t levelVertexCount[POP_MAXLEVEL + 1] = {0};
    size_t levelTriangleCount[POP_MAXLEVEL + 1] = {0};

    /* Every source vertex and face belongs to exactly one activation level.
     * Dense level vectors are substantially cheaper than the old
     * map<level, unordered_set<vertex>> representation on scan meshes: Lucy's
     * 14 million vertices otherwise spent more memory on hash nodes than on
     * the source coordinates themselves.  BoT indices are int-sized by
     * contract, so 32-bit remap entries are sufficient. */
    std::vector<uint32_t> triIndexMap;
    std::vector<uint8_t> vertexTriMinLevel;
    std::vector<std::vector<uint32_t>> levelTriVerts;
    std::vector<std::vector<uint32_t>> levelTris;

    size_t vertexCount = 0;
    const point_t *vertexArray = NULL;
    const vect_t *normalArray = NULL;
    size_t faceCount = 0;
    const int *faceArray = NULL;

    struct BObolMeshLodContext *context = NULL;
    struct bu_cache_txn *readTxn = NULL;
    struct bu_cache_txn *writeTxn = NULL;
    std::shared_lock<std::shared_mutex> readLock;
};

void
BObolPopState::triProcess(void)
{
    vertexTriMinLevel.assign(vertexCount,
	static_cast<uint8_t>(POP_MAXLEVEL - 1));
    levelTriVerts.resize(POP_MAXLEVEL);
    levelTris.resize(POP_MAXLEVEL);

    for (size_t faceIndex = 0; faceIndex < faceCount; faceIndex++) {
	BObolPopRec triangle[3];
	bool badFace = false;
	for (size_t cornerIndex = 0; cornerIndex < 3; cornerIndex++) {
	    int faceVertex = faceArray[3 * faceIndex + cornerIndex];
	    if (faceVertex < 0 || static_cast<size_t>(faceVertex) >= vertexCount) {
		bu_log("bad face %zd - skipping\n", faceIndex);
		badFace = true;
		break;
	    }
	    const auto quantize = [](fastf_t value, fastf_t minimum,
		    fastf_t maximum) -> unsigned short {
		const fastf_t extent = maximum - minimum;
		if (!(extent > 0.0))
		    return 0;
		const fastf_t normalized = (value - minimum) / extent;
		if (normalized <= 0.0)
		    return 0;
		if (normalized >= 1.0)
		    return USHRT_MAX;
		return static_cast<unsigned short>(
		    floor(normalized * static_cast<fastf_t>(USHRT_MAX)));
	    };
	    triangle[cornerIndex].x = quantize(
		vertexArray[faceVertex][X], minx, maxx);
	    triangle[cornerIndex].y = quantize(
		vertexArray[faceVertex][Y], miny, maxy);
	    triangle[cornerIndex].z = quantize(
		vertexArray[faceVertex][Z], minz, maxz);
	}
	if (badFace)
	    continue;

	/* A pair first separates at the level exposing its most-significant
	 * differing quantized bit.  A triangle activates once all three vertex
	 * pairs have separated.  Computing that level directly is equivalent to
	 * testing all 16 masks, but avoids up to 48 quotient comparisons for
	 * every face in very large meshes. */
	const auto pairLevel = [](const BObolPopRec &a,
		const BObolPopRec &b) -> int {
	    const unsigned int differing =
		static_cast<unsigned int>(a.x ^ b.x) |
		static_cast<unsigned int>(a.y ^ b.y) |
		static_cast<unsigned int>(a.z ^ b.z);
	    if (!differing)
		return POP_MAXLEVEL - 1;
#if defined(__GNUC__) || defined(__clang__)
	    const int mostSignificantBit =
		31 - __builtin_clz(differing);
#else
	    int mostSignificantBit = 0;
	    for (unsigned int bits = differing; bits >>= 1;)
		++mostSignificantBit;
#endif
	    return (POP_MAXLEVEL - 1) - mostSignificantBit;
	};
	const int level = std::max(pairLevel(triangle[0], triangle[1]),
	    std::max(pairLevel(triangle[1], triangle[2]),
		     pairLevel(triangle[0], triangle[2])));
	levelTris[level].push_back(static_cast<uint32_t>(faceIndex));

	for (size_t cornerIndex = 0; cornerIndex < 3; cornerIndex++) {
	    int faceVertex = faceArray[3 * faceIndex + cornerIndex];
	    if (vertexTriMinLevel[faceVertex] > level)
		vertexTriMinLevel[faceVertex] = static_cast<uint8_t>(level);
	}
    }

    for (size_t vertexIndex = 0; vertexIndex < vertexTriMinLevel.size();
	 vertexIndex++)
	levelTriVerts[vertexTriMinLevel[vertexIndex]].push_back(
	    static_cast<uint32_t>(vertexIndex));

    triIndexMap.resize(vertexCount);

    uint32_t outputVertexIndex = 0;
    for (const std::vector<uint32_t> &levelVertices : levelTriVerts) {
	for (uint32_t sourceVertexIndex : levelVertices) {
	    triIndexMap[sourceVertexIndex] = outputVertexIndex;
	    ++outputVertexIndex;
	}
    }

    int firstPopulatedLevel = -1;
    int lastPopulatedLevel = -1;
    for (size_t levelIndex = 0; levelIndex < levelTris.size(); levelIndex++) {
	if (levelTris[levelIndex].size()) {
	    if (firstPopulatedLevel < 0)
		firstPopulatedLevel = static_cast<int>(levelIndex);
	    lastPopulatedLevel = static_cast<int>(levelIndex);
	}
    }

    /* Face activation and coordinate precision are independent progress axes.
     * A mesh may have every face active at a low level while its vertices are
     * still snapped to a visibly coarse grid.  The last triangle-pop level is
     * therefore not a valid display terminal.  Persist the cumulative
     * topology through the full 16-bit precision range; levels with no new
     * faces or vertices are cheap and still re-snap the resident coordinates
     * to a finer grid. */
    minPopLevel = firstPopulatedLevel >= 0 ? firstPopulatedLevel : 0;
    maxPopLevel = lastPopulatedLevel >= 0 ? POP_MAXLEVEL - 1 : minPopLevel;
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

    if (!userKey) {
	struct bu_data_hash_state *state = bu_data_hash_create();
	static const char semantics[] = "BObol-PoP3D-format-10";
	bu_data_hash_update(state, semantics, sizeof(semantics));
	bu_data_hash_update(state, vertices, inputVertexCount * sizeof(point_t));
	bu_data_hash_update(state, faces, 3 * inputFaceCount * sizeof(int));
	const unsigned char normalFlag = normals ? 1u : 0u;
	bu_data_hash_update(state, &normalFlag, sizeof(normalFlag));
	if (normals)
	    bu_data_hash_update(state, normals,
		3 * inputFaceCount * sizeof(vect_t));
	const unsigned char cullFlag = cullBackfaces ? 1u : 0u;
	bu_data_hash_update(state, &cullFlag, sizeof(cullFlag));
	hash = bu_data_hash_val(state);
	bu_data_hash_destroy(state);
    } else {
	hash = userKey;
    }

    void *cacheData = NULL;
    size_t cacheSize = cacheGet(&cacheData, CACHE_POP_MAX_LEVEL);
    if (cacheSize && cacheData) {
	cacheDone();
	isValid = true;
	return;
    }
    cacheDone();

    currLevel = POP_MAXLEVEL - 1;
    for (int levelIndex = 0; levelIndex < POP_MAXLEVEL; levelIndex++)
	precomputedMasks.push_back(static_cast<unsigned short>(
				       pow(2, (POP_MAXLEVEL - levelIndex - 1))));

    vertexCount = inputVertexCount;
    vertexArray = vertices;
    normalArray = normals;
    faceCount = inputFaceCount;
    faceArray = faces;

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

    const int64_t classifyStarted = bu_gettime();
    triProcess();
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
				 unsigned long long key)
{
    context = ctx;
    if (!key)
	return;

    currLevel = -1;
    for (int levelIndex = 0; levelIndex < POP_MAXLEVEL; levelIndex++)
	precomputedMasks.push_back(static_cast<unsigned short>(
				       pow(2, (POP_MAXLEVEL - levelIndex - 1))));
    hash = key;

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
    if (!readMetadata(CACHE_POP_MAX_LEVEL, &maxPopLevel,
	    sizeof(maxPopLevel)) ||
	!readMetadata(CACHE_POP_MIN_LEVEL, &minPopLevel,
	    sizeof(minPopLevel)) ||
	!readMetadata(CACHE_SHADED_CULL_BACKFACES,
	    &shadedCullBackfaces, sizeof(shadedCullBackfaces)) ||
	!readMetadata(CACHE_VERTEX_COUNT, levelVertexCount,
	    sizeof(levelVertexCount)) ||
	!readMetadata(CACHE_TRI_COUNT, levelTriangleCount,
	    sizeof(levelTriangleCount))) {
	cacheDone();
	return;
    }

    {
	fastf_t minmax[6];
	const char *buffer = NULL;
	size_t bufferSize = cacheGet((void **)&buffer, CACHE_OBJ_BOUNDS);
	if (bufferSize != (sizeof(bbmin) + sizeof(bbmax) + sizeof(minmax))) {
	    cacheDone();
	    return;
	}
	memcpy(&bbmin, buffer, sizeof(bbmin));
	buffer += sizeof(bbmin);
	memcpy(&bbmax, buffer, sizeof(bbmax));
	buffer += sizeof(bbmax);
	memcpy(&minmax, buffer, sizeof(minmax));
	minx = minmax[0];
	miny = minmax[1];
	minz = minmax[2];
	maxx = minmax[3];
	maxy = minmax[4];
	maxz = minmax[5];
    }
    cacheDone();

    /* Metadata validation must not also read geometry.  The caller supplies
     * the view-selected target immediately after opening the handle, so
     * eagerly loading the minimum cut here doubled cache I/O for every
     * request above that cut. */
    isValid = true;
}

BObolPopState::~BObolPopState()
{
    if (readTxn || readLock.owns_lock())
	cacheDone();
    if (fullDetailFree) {
	fullDetailFree(detailData);
	detailData = NULL;
    } else if (fullDetailClear) {
	fullDetailClear(detailData);
	detailData = NULL;
    }
}

void
BObolPopState::updateSnappedPoints(int level)
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
	levelPoint(&snapped, &point, level);
	for (int coordIndex = 0; coordIndex < 3; coordIndex++)
	    lodTriPointsSnapped.push_back(snapped[coordIndex]);
    }
}

bool
BObolPopState::triPopLoad(int startLevel, int level,
			  bool materializeSnapped)
{
    struct bu_vls keyBuffer = BU_VLS_INIT_ZERO;

    /* A requested PoP prefix is one immutable cache snapshot.  Keep the LMDB
     * read transaction open while copying every vertex, triangle, and normal
     * chunk in the prefix.  Finishing a transaction after every component
     * made a 50k-asset warm scene perform hundreds of thousands of tiny read
     * transactions and reduced an otherwise memory-resident cache to roughly
     * a hundred meshes per second. */
    for (int levelIndex = startLevel + 1; levelIndex <= level; levelIndex++) {
	if (!levelVertexCount[levelIndex])
	    continue;
	bu_vls_sprintf(&keyBuffer, "%s%d", CACHE_VERT_LEVEL, levelIndex);
	fastf_t *buffer = NULL;
	size_t bufferSize = cacheGet((void **)&buffer, bu_vls_cstr(&keyBuffer));
	if (bufferSize != levelVertexCount[levelIndex] * sizeof(point_t)) {
	    cacheDone();
	    bu_vls_free(&keyBuffer);
	    return false;
	}
	lodTriPoints.insert(lodTriPoints.end(), &buffer[0],
			    &buffer[levelVertexCount[levelIndex] * 3]);
    }

    if (materializeSnapped)
	updateSnappedPoints(level);
    else
	lodTriPointsSnapped.clear();

    for (int levelIndex = startLevel + 1; levelIndex <= level; levelIndex++) {
	if (!levelTriangleCount[levelIndex])
	    continue;
	bu_vls_sprintf(&keyBuffer, "%s%d", CACHE_TRI_LEVEL, levelIndex);
	int *buffer = NULL;
	size_t bufferSize = cacheGet((void **)&buffer, bu_vls_cstr(&keyBuffer));
	if (bufferSize != levelTriangleCount[levelIndex] * 3 * sizeof(int)) {
	    cacheDone();
	    bu_vls_free(&keyBuffer);
	    return false;
	}
	lodTris.insert(lodTris.end(), &buffer[0],
		       &buffer[levelTriangleCount[levelIndex] * 3]);
    }

    for (int levelIndex = startLevel + 1; levelIndex <= level; levelIndex++) {
	if (!levelTriangleCount[levelIndex])
	    continue;
	bu_vls_sprintf(&keyBuffer, "%s%d", CACHE_VERTNORM_LEVEL, levelIndex);
	fastf_t *buffer = NULL;
	size_t bufferSize = cacheGet((void **)&buffer, bu_vls_cstr(&keyBuffer));
	if (bufferSize > 0 &&
	    bufferSize != levelTriangleCount[levelIndex] * sizeof(vect_t) * 3) {
	    cacheDone();
	    bu_vls_free(&keyBuffer);
	    return false;
	}
	if (bufferSize) {
	    lodTriNormals.insert(lodTriNormals.end(), &buffer[0],
				 &buffer[levelTriangleCount[levelIndex] * 3 * 3]);
	}
    }

    cacheDone();
    bu_vls_free(&keyBuffer);
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

void
BObolPopState::triPopTrim(int level, bool materializeSnapped)
{
    size_t keepVertexCount = 0;
    size_t keepFaceCount = 0;
    for (size_t levelIndex = 0; levelIndex <= static_cast<size_t>(level);
	 levelIndex++) {
	keepVertexCount += levelVertexCount[levelIndex];
	keepFaceCount += levelTriangleCount[levelIndex];
    }

    lodTriPoints.resize(keepVertexCount * 3);
    lodTriPoints.shrink_to_fit();
    lodTriNormals.resize(keepFaceCount * 3 * 3);
    lodTriNormals.shrink_to_fit();
    lodTris.resize(keepFaceCount * 3);
    lodTris.shrink_to_fit();

    if (materializeSnapped)
	updateSnappedPoints(level);
    else
	lodTriPointsSnapped.clear();
}

int
BObolPopState::getLevel(fastf_t viewLength)
{
    fastf_t delta = 0.01 * viewLength;
    point_t bmin;
    point_t bmax;
    VSET(bmin, minx, miny, minz);
    VSET(bmax, maxx, maxy, maxz);
    fastf_t bboxDiagonal = DIST_PNT_PNT(bmin, bmax);

    for (int levelIndex = 0; levelIndex < POP_MAXLEVEL; levelIndex++) {
	fastf_t diagSlice = bboxDiagonal / pow(2, levelIndex);
	if (diagSlice < delta)
	    return levelIndex;
    }
    return POP_MAXLEVEL - 1;
}

void
BObolPopState::hierarchyInfo(struct BObolMeshLodHierarchyInfo *info) const
{
    if (!info)
	return;
    info->min_level = minPopLevel;
    info->max_level = maxPopLevel;
    info->resident_level = currLevel;
    VSET(info->quantization_min, minx, miny, minz);
    VSET(info->quantization_max, maxx, maxy, maxz);
    size_t points = 0;
    size_t faces = 0;
    for (int level = 0; level < BOBOL_MESH_LOD_LEVEL_COUNT; ++level) {
	points += levelVertexCount[level];
	faces += levelTriangleCount[level];
	info->point_count[level] = points;
	info->face_count[level] = faces;
    }
}

bool
BObolPopState::setupFullDetail(void)
{
    if (!lod || !fullDetailSetup)
	return false;

    struct BObolMeshLodDetail detail;
    bobol_mesh_lod_detail_init(&detail);
    if (fullDetailSetup(&detail, detailData) != 0)
	return false;

    if (!mesh_lod_arrays_validate(detail.faces, detail.face_count,
				  detail.points, detail.point_count,
				  detail.points_orig, detail.point_orig_count,
				  NULL, NULL, NULL)) {
	if (fullDetailClear)
	    fullDetailClear(detailData);
	return false;
    }

    const size_t indexCount = detail.face_count * 3;
    if ((detail.normals && detail.normal_count != indexCount) ||
	(!detail.normals && detail.normal_count != 0)) {
	if (fullDetailClear)
	    fullDetailClear(detailData);
	return false;
    }

    lod->faces = detail.faces;
    lod->fcnt = static_cast<int>(detail.face_count);
    lod->points = detail.points;
    lod->pcnt = static_cast<int>(detail.point_count);
    lod->points_orig = detail.points_orig;
    lod->porig_cnt = static_cast<int>(detail.point_orig_count);
    lod->normals = detail.normals;
    return true;
}

void
BObolPopState::clearFullDetail(void)
{
    if (fullDetailClear)
	fullDetailClear(detailData);
    mesh_lod_active_data_clear(lod);
}

bool
BObolPopState::setLevel(int level, bool materializeSnapped)
{
    if (level > maxPopLevel && !fullDetailSetup)
	level = maxPopLevel;

    if (level == currLevel && !forceUpdate) {
	if (materializeSnapped && lodTriPointsSnapped.empty() &&
	    !lodTriPoints.empty())
	    updateSnappedPoints(level);
	return true;
    }

    if (forceUpdate) {
	forceUpdate = false;
	clearFullDetail();
	shrinkMemory();
	currLevel = -1;
	return setLevel(level, materializeSnapped);
    }

    if (level > currLevel && level <= maxPopLevel) {
	if (!lodTriPoints.size()) {
	    if (!triPopLoad(-1, level, materializeSnapped))
		return false;
	} else if (!triPopLoad(currLevel, level, materializeSnapped)) {
	    return false;
	}
	/* A retained exact-only load may be followed by a legacy snapped load at
	 * the same or a higher level.  Populate that derived array on demand
	 * without reopening any cache chunks. */
	if (materializeSnapped && lodTriPointsSnapped.empty())
	    updateSnappedPoints(level);
    }

    if (level < currLevel && level <= maxPopLevel &&
	currLevel <= maxPopLevel) {
	if (!lodTriPoints.size()) {
	    if (!triPopLoad(-1, level, materializeSnapped))
		return false;
	} else {
	    triPopTrim(level, materializeSnapped);
	}
    }

    if (level < currLevel && level <= maxPopLevel &&
	currLevel > maxPopLevel) {
	clearFullDetail();
	if (!triPopLoad(-1, level, materializeSnapped))
	    return false;
    }

    if (level > currLevel && level > maxPopLevel &&
	currLevel <= maxPopLevel) {
	lodTriPointsSnapped.clear();
	lodTriPointsSnapped.shrink_to_fit();
	lodTriPoints.clear();
	lodTriPoints.shrink_to_fit();
	lodTriNormals.clear();
	lodTriNormals.shrink_to_fit();
	lodTris.clear();
	lodTris.shrink_to_fit();

	if (!setupFullDetail()) {
	    currLevel = -1;
	    return false;
	}
    }

    currLevel = level;
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
    std::string keystr = std::to_string(hash) + std::string(":") +
			 std::string(component);
    size_t wsize = bu_cache_write(const_cast<void *>(data), size,
				  keystr.c_str(), context->i->lodCache,
				  &writeTxn);
    return (wsize > 0);
}

size_t
BObolPopState::cacheGet(void **data, const char *component)
{
    if (context && context->i && context->i->accessMutex &&
	!readLock.owns_lock())
	readLock =
	    std::shared_lock<std::shared_mutex>(*context->i->accessMutex);
    std::string keystr = std::to_string(hash) + std::string(":") +
			 std::string(component);
    return bu_cache_get(data, keystr.c_str(), context->i->lodCache,
	&readTxn);
}

void
BObolPopState::cacheDone(void)
{
    bu_cache_get_done(&readTxn);
    if (readLock.owns_lock())
	readLock.unlock();
}

bool
BObolPopState::cacheTri(void)
{
    {
	std::stringstream stream;
	stream.write(reinterpret_cast<const char *>(&maxPopLevel),
		     sizeof(maxPopLevel));
	if (!cacheWrite(CACHE_POP_MAX_LEVEL, stream))
	    return false;
    }

    {
	std::stringstream stream;
	stream.write(reinterpret_cast<const char *>(&minPopLevel),
		     sizeof(minPopLevel));
	if (!cacheWrite(CACHE_POP_MIN_LEVEL, stream))
	    return false;
    }

    {
	std::stringstream stream;
	stream.write(reinterpret_cast<const char *>(&shadedCullBackfaces),
		     sizeof(shadedCullBackfaces));
	if (!cacheWrite(CACHE_SHADED_CULL_BACKFACES, stream))
	    return false;
    }

    {
	std::stringstream stream;
	for (size_t levelIndex = 0; levelIndex <= POP_MAXLEVEL; levelIndex++) {
	    size_t count = 0;
	    if (levelIndex >= levelTriVerts.size()) {
		stream.write(reinterpret_cast<const char *>(&count),
			     sizeof(count));
		continue;
	    }
	    if (static_cast<int>(levelIndex) > maxPopLevel ||
		levelTriVerts[levelIndex].empty()) {
		stream.write(reinterpret_cast<const char *>(&count),
			     sizeof(count));
		continue;
	    }
	    count = levelTriVerts[levelIndex].size();
	    stream.write(reinterpret_cast<const char *>(&count), sizeof(count));
	}
	if (!cacheWrite(CACHE_VERTEX_COUNT, stream))
	    return false;
    }

    {
	std::stringstream stream;
	for (size_t levelIndex = 0; levelIndex <= POP_MAXLEVEL; levelIndex++) {
	    size_t count = 0;
	    if (static_cast<int>(levelIndex) > maxPopLevel ||
		!levelTris[levelIndex].size()) {
		stream.write(reinterpret_cast<const char *>(&count),
			     sizeof(count));
		continue;
	    }
	    count = levelTris[levelIndex].size();
	    stream.write(reinterpret_cast<const char *>(&count), sizeof(count));
	}
	if (!cacheWrite(CACHE_TRI_COUNT, stream))
	    return false;
    }

    struct bu_vls keyBuffer = BU_VLS_INIT_ZERO;

    for (int levelIndex = 0; levelIndex <= maxPopLevel;
	 levelIndex++) {
	if (levelIndex >= static_cast<int>(levelTriVerts.size()) ||
	    levelTriVerts[levelIndex].empty())
	    continue;
	std::vector<fastf_t> points;
	points.resize(levelTriVerts[levelIndex].size() * 3);
	size_t outputIndex = 0;
	for (uint32_t sourceVertexIndex : levelTriVerts[levelIndex]) {
	    points[outputIndex++] = vertexArray[sourceVertexIndex][X];
	    points[outputIndex++] = vertexArray[sourceVertexIndex][Y];
	    points[outputIndex++] = vertexArray[sourceVertexIndex][Z];
	}
	bu_vls_sprintf(&keyBuffer, "%s%d", CACHE_VERT_LEVEL, levelIndex);
	if (!cacheWriteData(bu_vls_cstr(&keyBuffer), points.data(),
		points.size() * sizeof(fastf_t))) {
	    bu_vls_free(&keyBuffer);
	    return false;
	}
    }

    for (int levelIndex = 0; levelIndex <= maxPopLevel;
	 levelIndex++) {
	if (levelTris[levelIndex].empty())
	    continue;
	std::vector<int> triangles;
	triangles.resize(levelTris[levelIndex].size() * 3);
	size_t outputIndex = 0;
	for (uint32_t sourceFaceIndex : levelTris[levelIndex]) {
	    triangles[outputIndex++] = static_cast<int>(
		triIndexMap[faceArray[3 * sourceFaceIndex + 0]]);
	    triangles[outputIndex++] = static_cast<int>(
		triIndexMap[faceArray[3 * sourceFaceIndex + 1]]);
	    triangles[outputIndex++] = static_cast<int>(
		triIndexMap[faceArray[3 * sourceFaceIndex + 2]]);
	}
	bu_vls_sprintf(&keyBuffer, "%s%d", CACHE_TRI_LEVEL, levelIndex);
	if (!cacheWriteData(bu_vls_cstr(&keyBuffer), triangles.data(),
		triangles.size() * sizeof(int))) {
	    bu_vls_free(&keyBuffer);
	    return false;
	}
    }

    if (normalArray) {
	for (int levelIndex = 0; levelIndex <= maxPopLevel;
	     levelIndex++) {
	    if (levelTris[levelIndex].empty())
		continue;
	    std::vector<fastf_t> normals;
	    normals.resize(levelTris[levelIndex].size() * 9);
	    size_t outputIndex = 0;
	    for (uint32_t sourceFaceIndex : levelTris[levelIndex]) {
		for (int cornerIndex = 0; cornerIndex < 3; cornerIndex++) {
		    const size_t normalIndex =
			3 * static_cast<size_t>(sourceFaceIndex) + cornerIndex;
		    normals[outputIndex++] = normalArray[normalIndex][X];
		    normals[outputIndex++] = normalArray[normalIndex][Y];
		    normals[outputIndex++] = normalArray[normalIndex][Z];
		}
	    }
	    bu_vls_sprintf(&keyBuffer, "%s%d", CACHE_VERTNORM_LEVEL,
			   levelIndex);
	    if (!cacheWriteData(bu_vls_cstr(&keyBuffer), normals.data(),
		    normals.size() * sizeof(fastf_t))) {
		bu_vls_free(&keyBuffer);
		return false;
	    }
	}
    }

    bu_vls_free(&keyBuffer);
    return true;
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
	std::stringstream stream;
	stream.write(reinterpret_cast<const char *>(&bbmin), sizeof(bbmin));
	stream.write(reinterpret_cast<const char *>(&bbmax), sizeof(bbmax));
	stream.write(reinterpret_cast<const char *>(&minx), sizeof(minx));
	stream.write(reinterpret_cast<const char *>(&miny), sizeof(miny));
	stream.write(reinterpret_cast<const char *>(&minz), sizeof(minz));
	stream.write(reinterpret_cast<const char *>(&maxx), sizeof(maxx));
	stream.write(reinterpret_cast<const char *>(&maxy), sizeof(maxy));
	stream.write(reinterpret_cast<const char *>(&maxz), sizeof(maxz));
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

int
BObolPopState::toLevel(int value, int level)
{
    return static_cast<int>(floor(value / double(precomputedMasks[level])));
}

fastf_t
BObolPopState::snap(fastf_t value, fastf_t min, fastf_t max, int level)
{
    if (!(max > min))
	return value;
    unsigned int vf = static_cast<unsigned int>(
			  floor((value - min) / (max - min) * USHRT_MAX));
    int lv = static_cast<int>(floor(vf / double(precomputedMasks[level])));
    unsigned int vc = static_cast<unsigned int>(
			  ceil((value - min) / (max - min) * USHRT_MAX));
    int hc = static_cast<int>(ceil(vc / double(precomputedMasks[level])));
    fastf_t snapped =
	(static_cast<fastf_t>(lv) + static_cast<fastf_t>(hc)) * 0.5 *
	double(precomputedMasks[level]);
    return ((snapped / USHRT_MAX) * (max - min)) + min;
}

void
BObolPopState::levelPoint(point_t *out, const point_t *point, int level)
{
    /* The terminal cut is the stable exact state.  All source vertices and
     * faces are already present in the cumulative cache at this point, so
     * retaining the original coordinate costs no additional I/O and avoids
     * imposing a permanent 16-bit quantization error on very large extents. */
    if (level >= POP_MAXLEVEL - 1) {
	VMOVE(*out, *point);
	return;
    }
    fastf_t nx = snap((*point)[X], minx, maxx, level);
    fastf_t ny = snap((*point)[Y], miny, maxy, level);
    fastf_t nz = snap((*point)[Z], minz, maxz, level);
    VSET(*out, nx, ny, nz);
}

bool
BObolPopState::isEqual(BObolPopRec r1, BObolPopRec r2, int level)
{
    bool equalX = (toLevel(r1.x, level) == toLevel(r2.x, level));
    bool equalY = (toLevel(r1.y, level) == toLevel(r2.y, level));
    bool equalZ = (toLevel(r1.z, level) == toLevel(r2.z, level));
    return (equalX && equalY && equalZ);
}

bool
BObolPopState::triDegenerate(BObolPopRec r0,
			       BObolPopRec r1,
			       BObolPopRec r2,
			       int level)
{
    return isEqual(r0, r1, level) || isEqual(r1, r2, level) ||
	   isEqual(r0, r2, level);
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
		unsigned long long key)
{
    if (!context || !key)
	return NULL;

    BObolPopState *state = new (std::nothrow) BObolPopState(context, key);
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

static int
mesh_lod_level(struct BObolMeshLod *lod, int level, int reset,
	       bool materializeSnapped = true)
{
    if (!lod || !lod->state)
	return -1;
    BObolPopState *state = lod->state;
    if (level < 0)
	return state->currLevel;

    if (reset)
	state->forceUpdate = true;
    if (!state->setLevel(level, materializeSnapped)) {
	mesh_lod_active_data_clear(lod);
	return -1;
    }

    if (state->currLevel <= state->maxPopLevel)
	mesh_lod_active_pop_data_publish(lod, state);

    return state->currLevel;
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

    if (context && context->i && context->i->accessMutex)
	lock =
	    std::unique_lock<std::shared_mutex>(*context->i->accessMutex);

    if (context && key) {
	mesh_lod_cache_del(context, key, CACHE_POP_MAX_LEVEL);
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

static int
mesh_lod_bot_detail_setup(struct BObolMeshLodDetail *detail, void *cbData)
{
    if (!detail || !cbData)
	return -1;

    struct BObolMeshLodBotDetail *callbackData =
	    static_cast<struct BObolMeshLodBotDetail *>(cbData);
    if (!callbackData->dbip || !callbackData->dp)
	return -1;

    if (callbackData->intern)
	return 0;

    BU_GET(callbackData->intern, struct rt_db_internal);
    RT_DB_INTERNAL_INIT(callbackData->intern);
    if (rt_db_get_internal(callbackData->intern, callbackData->dp,
			   callbackData->dbip, NULL) < 0) {
	BU_PUT(callbackData->intern, struct rt_db_internal);
	callbackData->intern = NULL;
	return -1;
    }

    if (callbackData->intern->idb_minor_type != DB5_MINORTYPE_BRLCAD_BOT) {
	rt_db_free_internal(callbackData->intern);
	BU_PUT(callbackData->intern, struct rt_db_internal);
	callbackData->intern = NULL;
	return -1;
    }

    struct rt_bot_internal *bot =
	    static_cast<struct rt_bot_internal *>(callbackData->intern->idb_ptr);
    RT_BOT_CK_MAGIC(bot);

    int faceCount = 0;
    int pointCount = 0;
    int pointOrigCount = 0;
    if (!mesh_lod_arrays_validate(bot->faces, bot->num_faces,
				  reinterpret_cast<const point_t *>(bot->vertices),
				  bot->num_vertices,
				  reinterpret_cast<const point_t *>(bot->vertices),
				  bot->num_vertices,
				  &faceCount,
				  &pointCount,
				  &pointOrigCount)) {
	rt_db_free_internal(callbackData->intern);
	BU_PUT(callbackData->intern, struct rt_db_internal);
	callbackData->intern = NULL;
	return -1;
    }

    detail->faces = bot->faces;
    detail->face_count = static_cast<size_t>(faceCount);
    detail->points = reinterpret_cast<const point_t *>(bot->vertices);
    detail->point_count = static_cast<size_t>(pointCount);
    detail->points_orig = reinterpret_cast<const point_t *>(bot->vertices);
    detail->point_orig_count = static_cast<size_t>(pointOrigCount);
    detail->normals = NULL;
    detail->normal_count = 0;

    return 0;
}

static int
mesh_lod_bot_detail_clear(void *cbData)
{
    struct BObolMeshLodBotDetail *callbackData =
	    static_cast<struct BObolMeshLodBotDetail *>(cbData);

    if (callbackData && callbackData->intern) {
	rt_db_free_internal(callbackData->intern);
	BU_PUT(callbackData->intern, struct rt_db_internal);
	callbackData->intern = NULL;
    }

    return 0;
}

static int
mesh_lod_bot_detail_free(void *cbData)
{
    mesh_lod_bot_detail_clear(cbData);
    if (cbData) {
	struct BObolMeshLodBotDetail *callbackData =
		static_cast<struct BObolMeshLodBotDetail *>(cbData);
	BU_PUT(callbackData, struct BObolMeshLodBotDetail);
    }

    return 0;
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

int
bobol_mesh_lod_cache_refresh(struct db_i *dbip,
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
    if (rt_db_get_internal(&dbintern, dp, dbip, NULL) < 0) {
	mesh_lod_context_destroy(context);
	return BRLCAD_ERROR;
    }
    if (dbintern.idb_minor_type != DB5_MINORTYPE_BRLCAD_BOT) {
	rt_db_free_internal(&dbintern);
	mesh_lod_context_destroy(context);
	return BRLCAD_ERROR;
    }

    struct rt_bot_internal *bot =
	    static_cast<struct rt_bot_internal *>(dbintern.idb_ptr);
    RT_BOT_CK_MAGIC(bot);

    if (!mesh_lod_arrays_validate(bot->faces, bot->num_faces,
				  reinterpret_cast<const point_t *>(bot->vertices),
				  bot->num_vertices,
				  reinterpret_cast<const point_t *>(bot->vertices),
				  bot->num_vertices, NULL, NULL, NULL)) {
	rt_db_free_internal(&dbintern);
	if (status)
	    *status = current;
	mesh_lod_context_destroy(context);
	return BRLCAD_ERROR;
    }

    /* Obol's cull-safe contract is always exterior CCW.  Begin with the
     * declared winding, then validate actual topology.  Imported meshes such
     * as Stanford Lucy are often marked unoriented despite having one closed,
     * consistently wound shell; normalize those from their signed volume
     * instead of permanently forcing two-sided shaded rasterization. */
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

    bool cullBackfaces = false;
    if (bot->mode == RT_BOT_SOLID &&
	bot->num_vertices <= static_cast<size_t>(INT_MAX) &&
	bot->num_faces <= static_cast<size_t>(INT_MAX) &&
	bg_trimesh_solid2(static_cast<int>(bot->num_vertices),
	    static_cast<int>(bot->num_faces), bot->vertices,
	    const_cast<int *>(cacheFaces), NULL) == 0) {
	const bool oneComponent = mesh_lod_faces_one_component(
	    cacheFaces, bot->num_faces, bot->num_vertices);
	if (oneComponent) {
	    const int windingSign = mesh_lod_closed_winding_sign(
		reinterpret_cast<const point_t *>(bot->vertices),
		bot->num_vertices, cacheFaces, bot->num_faces);
	    if (windingSign != 0) {
		if (windingSign < 0) {
		    flippedWinding = !flippedWinding;
		    applyWinding();
		}
		cullBackfaces = true;
	    }
	} else if (bot->orientation != RT_BOT_UNORIENTED) {
	    /* Multiple shells can encode cavities, so their signed volumes
	     * cannot be normalized independently.  An explicit orientation is
	     * the only available material-side contract in that case. */
	    cullBackfaces = true;
	}
    }

    std::vector<fastf_t> authoredNormals;
    const point_t *botVertices =
	reinterpret_cast<const point_t *>(bot->vertices);
    const vect_t *botNormals =
	mesh_lod_bot_authored_normals(
	    authoredNormals, bot, cacheFaces, flippedWinding) ?
	reinterpret_cast<const vect_t *>(authoredNormals.data()) : NULL;

    unsigned long long key = mesh_lod_cache_generate(
				 context, botVertices, bot->num_vertices,
				 botNormals, cacheFaces, bot->num_faces, 0,
				 cullBackfaces);
    if (!key || mesh_lod_key_put(context, dp->d_namep, key) != 0) {
	rt_db_free_internal(&dbintern);
	mesh_lod_context_destroy(context);
	return BRLCAD_ERROR;
    }

    rt_db_free_internal(&dbintern);

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

    if (current.has_cache_key &&
	(current.cache_key != userKey || !current.has_cached_payload)) {
	current.cleared_cache_entry = 1;
	current.cleared_cache_key = current.cache_key;
	mesh_lod_cache_clear_context(context, current.cache_key);
	current.cache_key = 0;
	current.has_cache_key = 0;
	current.has_cached_payload = 0;
	current.stale_cache_entry = 0;
    }

    if (current.has_cache_key && current.has_cached_payload) {
	if (status)
	    *status = current;
	mesh_lod_context_destroy(context);
	return BRLCAD_OK;
    }

    const bool cullBackfaces = shadedCullBackfaces &&
	vertexCount <= static_cast<size_t>(INT_MAX) &&
	faceCount <= static_cast<size_t>(INT_MAX) &&
	bg_trimesh_solid2(static_cast<int>(vertexCount),
	    static_cast<int>(faceCount),
	    const_cast<fastf_t *>(reinterpret_cast<const fastf_t *>(vertices)),
	    const_cast<int *>(faces), NULL) == 0;

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

    if (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) {
	struct BObolMeshLodBotDetail *callbackData;
	BU_GET(callbackData, struct BObolMeshLodBotDetail);
	memset(callbackData, 0, sizeof(*callbackData));
	callbackData->dbip = dbip;
	callbackData->dp = dp;
	if (!bobol_mesh_lod_detail_callbacks_set(
		lod, mesh_lod_bot_detail_setup, mesh_lod_bot_detail_clear,
		mesh_lod_bot_detail_free, callbackData))
	    mesh_lod_bot_detail_free(callbackData);
    }

    return lod;
}

unsigned long long
bobol_mesh_lod_cache_key_get(const struct BObolMeshLod *lod)
{
    return lod && lod->state ? lod->state->hash : 0;
}

int
bobol_mesh_lod_load_level(struct BObolMeshLod *lod, int level, int reset)
{
    return mesh_lod_level(lod, level, reset);
}

int
bobol_mesh_lod_load_display_level(struct BObolMeshLod *lod, int level,
	int reset)
{
    if (!lod || !lod->state)
	return -1;
    if (level < 0)
	level = lod->state->minPopLevel;
    level = std::max(level, lod->state->minPopLevel);
    level = std::min(level, lod->state->maxPopLevel);
    return mesh_lod_level(lod, level, reset);
}

int
bobol_mesh_lod_load_resident_level(struct BObolMeshLod *lod, int level,
	int reset)
{
    if (!lod || !lod->state)
	return -1;
    if (level < 0)
	level = lod->state->minPopLevel;
    level = std::max(level, lod->state->minPopLevel);
    level = std::min(level, lod->state->maxPopLevel);
    return mesh_lod_level(lod, level, reset, false);
}

int
bobol_mesh_lod_load_view(struct BObolMeshLod *lod,
			   const struct bv_view_info *info,
			   int reset)
{
    if (!lod || !lod->state)
	return -1;

    struct bv_view_info sanitized = BV_VIEW_INFO_INIT;
    if (info)
	sanitized = *info;
    bv_view_info_sanitize(&sanitized);

    fastf_t policyScale = sanitized.lod.scale;
    if (sanitized.size <= SMALL_FASTF)
	sanitized.size = 1.0;
    if (policyScale <= SMALL_FASTF)
	policyScale = 1.0;

    BObolPopState *state = lod->state;
    int viewLevel = static_cast<int>(
			static_cast<double>(state->getLevel(sanitized.size)) * policyScale);
    viewLevel = (viewLevel < 0) ? 0 : viewLevel;
    viewLevel = (viewLevel >= POP_MAXLEVEL) ? POP_MAXLEVEL - 1 : viewLevel;
    /* A view/display request loads only a cumulative PoP cut.  The complete
     * source arrays remain behind the explicit exact/full-detail API. */
    viewLevel = std::max(viewLevel, state->minPopLevel);
    viewLevel = std::min(viewLevel, state->maxPopLevel);

    return mesh_lod_level(lod, viewLevel, reset);
}

int
bobol_mesh_lod_current_level(const struct BObolMeshLod *lod)
{
    return (lod && lod->state) ? lod->state->currLevel : -1;
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

    info->active_level = lod->state->currLevel;
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
	 lod->state->currLevel <= lod->state->maxPopLevel) ? 1 : 0;
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
    return info->min_level >= 0 && info->max_level >= info->min_level;
}

void
bobol_mesh_lod_detail_init(struct BObolMeshLodDetail *detail)
{
    if (detail)
	memset(detail, 0, sizeof(*detail));
}

int
bobol_mesh_lod_detail_callbacks_set(
    struct BObolMeshLod *lod,
    BObolMeshLodDetailSetupCallback setupCallback,
    BObolMeshLodDetailClearCallback clearCallback,
    BObolMeshLodDetailFreeCallback freeCallback,
    void *cbData)
{
    if (!lod || !lod->state || !setupCallback)
	return 0;

    BObolPopState *state = lod->state;
    const int hadCallbacks =
	(state->fullDetailSetup || state->fullDetailClear ||
	 state->fullDetailFree) ? 1 : 0;
    if (state->fullDetailFree)
	state->fullDetailFree(state->detailData);
    else if (state->fullDetailClear)
	state->fullDetailClear(state->detailData);

    if (hadCallbacks) {
	if (state->currLevel < 0) {
	    mesh_lod_active_data_clear(lod);
	} else if (state->currLevel > state->maxPopLevel) {
	    state->currLevel = -1;
	    mesh_lod_active_data_clear(lod);
	} else {
	    mesh_lod_active_pop_data_publish(lod, state);
	}
    }

    state->fullDetailSetup = setupCallback;
    state->fullDetailClear = clearCallback;
    state->fullDetailFree = freeCallback;
    state->detailData = cbData;
    return 1;
}

void
bobol_mesh_lod_detail_callbacks_clear(struct BObolMeshLod *lod)
{
    if (!lod || !lod->state)
	return;

    BObolPopState *state = lod->state;
    if (state->fullDetailFree)
	state->fullDetailFree(state->detailData);
    else if (state->fullDetailClear)
	state->fullDetailClear(state->detailData);
    state->fullDetailSetup = NULL;
    state->fullDetailClear = NULL;
    state->fullDetailFree = NULL;
    state->detailData = NULL;
    if (state->currLevel > state->maxPopLevel) {
	state->currLevel = -1;
	mesh_lod_active_data_clear(lod);
    }
}

void
bobol_mesh_lod_memshrink(struct BObolMeshLod *lod)
{
    if (!lod || !lod->state)
	return;

    if (lod->state->currLevel > lod->state->maxPopLevel)
	return;

    lod->state->shrinkMemory();
    lod->state->forceUpdate = true;
    mesh_lod_active_data_clear(lod);
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
