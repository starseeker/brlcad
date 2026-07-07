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

#include "brlobol/draw_cache.h"
#include "brlobol/mesh_lod_cache.h"

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
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#define POP_MAXLEVEL 16
#define POP_CACHEDIR BRLOBOL_DRAW_CACHE_DIR
#define MBUMP 1.01
#define CACHE_CURRENT_FORMAT 4

#define CACHE_POP_MAX_LEVEL "th"
#define CACHE_POP_SWITCH_LEVEL "sw"
#define CACHE_VERTEX_COUNT "vc"
#define CACHE_TRI_COUNT "tc"
#define CACHE_OBJ_BOUNDS "bb"
#define CACHE_VERT_LEVEL "v"
#define CACHE_VERTNORM_LEVEL "vn"
#define CACHE_TRI_LEVEL "t"

struct BRLObolMeshLodContextInternal {
    struct bu_cache *lodCache;
    struct bu_cache *nameCache;
    struct bu_vls *fname;
    char *registryKey;
    std::mutex *accessMutex;
};

struct BRLObolMeshLodContext {
    struct BRLObolMeshLodContextInternal *i;
    size_t refs;
};

class BRLObolPopState;

struct BRLObolMeshLod {
    struct BRLObolMeshLodContext *context;
    BRLObolPopState *state;

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

struct BRLObolMeshLodBotDetail {
    struct db_i *dbip;
    struct directory *dp;
    struct rt_db_internal *intern;
};

class BRLObolPopRec
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
	sem = bu_semaphore_register("BRLOBOL_MESH_LOD_CACHE_WRITE");
    return sem;
}

static std::mutex &
mesh_lod_context_registry_mutex(void)
{
    static std::mutex m;
    return m;
}

static std::map<std::string, struct BRLObolMeshLodContext *> &
mesh_lod_context_registry(void)
{
    static std::map<std::string, struct BRLObolMeshLodContext *> registry;
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

static void
mesh_lod_generate_crease_normals(std::vector<fastf_t> &normalStorage,
				 const point_t *points,
				 size_t pointCount,
				 const int *faces,
				 size_t faceCount)
{
    normalStorage.clear();
    if (!points || !pointCount || !faces || !faceCount)
	return;

    std::vector<fastf_t> faceNormals(faceCount * 3, 0.0);
    std::vector<size_t> vertexFaceCounts(pointCount, 0);

    for (size_t faceIndex = 0; faceIndex < faceCount; faceIndex++) {
	const int *face = &faces[3 * faceIndex];
	const int ia = face[0];
	const int ib = face[1];
	const int ic = face[2];
	if (ia < 0 || ib < 0 || ic < 0 ||
	    static_cast<size_t>(ia) >= pointCount ||
	    static_cast<size_t>(ib) >= pointCount ||
	    static_cast<size_t>(ic) >= pointCount)
	    return;

	vect_t ab;
	vect_t ac;
	vect_t normal;
	VSUB2(ab, points[ib], points[ia]);
	VSUB2(ac, points[ic], points[ia]);
	VCROSS(normal, ab, ac);
	if (MAGNITUDE(normal) > SMALL_FASTF)
	    VUNITIZE(normal);
	else
	    VSET(normal, 0.0, 0.0, 1.0);

	faceNormals[3 * faceIndex + X] = normal[X];
	faceNormals[3 * faceIndex + Y] = normal[Y];
	faceNormals[3 * faceIndex + Z] = normal[Z];
	vertexFaceCounts[static_cast<size_t>(ia)]++;
	vertexFaceCounts[static_cast<size_t>(ib)]++;
	vertexFaceCounts[static_cast<size_t>(ic)]++;
    }

    std::vector<size_t> vertexFaceOffsets(pointCount + 1, 0);
    for (size_t vertexIndex = 0; vertexIndex < pointCount; vertexIndex++)
	vertexFaceOffsets[vertexIndex + 1] =
	    vertexFaceOffsets[vertexIndex] + vertexFaceCounts[vertexIndex];

    std::vector<size_t> vertexFaceCursor = vertexFaceOffsets;
    std::vector<size_t> vertexFaces(faceCount * 3, 0);
    for (size_t faceIndex = 0; faceIndex < faceCount; faceIndex++) {
	const int *face = &faces[3 * faceIndex];
	for (int cornerIndex = 0; cornerIndex < 3; cornerIndex++) {
	    const size_t vertexIndex = static_cast<size_t>(face[cornerIndex]);
	    vertexFaces[vertexFaceCursor[vertexIndex]++] = faceIndex;
	}
    }

    normalStorage.assign(faceCount * 9, 0.0);
    const fastf_t creaseCos = 0.5;
    for (size_t faceIndex = 0; faceIndex < faceCount; faceIndex++) {
	const fastf_t *baseNormal = &faceNormals[3 * faceIndex];
	for (int cornerIndex = 0; cornerIndex < 3; cornerIndex++) {
	    const size_t vertexIndex =
		static_cast<size_t>(faces[3 * faceIndex + cornerIndex]);
	    vect_t smoothNormal;
	    VSET(smoothNormal, 0.0, 0.0, 0.0);
	    for (size_t i = vertexFaceOffsets[vertexIndex];
		 i < vertexFaceOffsets[vertexIndex + 1]; i++) {
		const size_t adjacentFace = vertexFaces[i];
		const fastf_t *adjacentNormal =
		    &faceNormals[3 * adjacentFace];
		const fastf_t dot =
		    baseNormal[X] * adjacentNormal[X] +
		    baseNormal[Y] * adjacentNormal[Y] +
		    baseNormal[Z] * adjacentNormal[Z];
		if (dot >= creaseCos) {
		    smoothNormal[X] += adjacentNormal[X];
		    smoothNormal[Y] += adjacentNormal[Y];
		    smoothNormal[Z] += adjacentNormal[Z];
		}
	    }
	    if (MAGNITUDE(smoothNormal) > SMALL_FASTF)
		VUNITIZE(smoothNormal);
	    else
		VMOVE(smoothNormal, baseNormal);

	    const size_t normalIndex = 3 * (3 * faceIndex + cornerIndex);
	    normalStorage[normalIndex + X] = smoothNormal[X];
	    normalStorage[normalIndex + Y] = smoothNormal[Y];
	    normalStorage[normalIndex + Z] = smoothNormal[Z];
	}
    }
}

static void
mesh_lod_active_data_clear(struct BRLObolMeshLod *lod)
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
mesh_lod_context_close(struct BRLObolMeshLodContext *context)
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
	delete context->i->accessMutex;
	BU_PUT(context->i, struct BRLObolMeshLodContextInternal);
    }
    BU_PUT(context, struct BRLObolMeshLodContext);
}

static void
mesh_lod_context_destroy(struct BRLObolMeshLodContext *context)
{
    if (!context)
	return;

    {
	std::lock_guard<std::mutex> guard(mesh_lod_context_registry_mutex());
	if (context->refs > 1) {
	    context->refs--;
	    return;
	}
	if (context->i && context->i->registryKey)
	    mesh_lod_context_registry().erase(context->i->registryKey);
	context->refs = 0;
    }

    mesh_lod_context_close(context);
}

static void
mesh_lod_cache_clear_context(struct BRLObolMeshLodContext *context,
			     unsigned long long key);

static struct BRLObolMeshLodContext *
mesh_lod_context_create(const char *name)
{
    if (!name)
	return NULL;

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s", bu_path_normalize(name));
    if (bu_vls_strlen(&fname) < 10)
	bu_vls_printf(&fname, "GGGGGGGGGGGGG");
    unsigned long long hash = bu_data_hash(bu_vls_cstr(&fname),
					   bu_vls_strlen(&fname) * sizeof(char));
    bu_path_component(&fname, bu_path_normalize(name), BU_PATH_BASENAME_EXTLESS);
    bu_vls_printf(&fname, "_%llu", hash);

    std::string registryKey(bu_vls_cstr(&fname));
    std::lock_guard<std::mutex> guard(mesh_lod_context_registry_mutex());
    {
	auto it = mesh_lod_context_registry().find(registryKey);
	if (it != mesh_lod_context_registry().end()) {
	    it->second->refs++;
	    bu_vls_free(&fname);
	    return it->second;
	}
    }

    struct BRLObolMeshLodContext *context;
    BU_GET(context, struct BRLObolMeshLodContext);
    BU_GET(context->i, struct BRLObolMeshLodContextInternal);
    context->refs = 1;
    struct BRLObolMeshLodContextInternal *internal = context->i;
    BU_GET(internal->fname, struct bu_vls);
    bu_vls_init(internal->fname);
    bu_vls_sprintf(internal->fname, "%s", bu_vls_cstr(&fname));
    internal->registryKey = bu_strdup(registryKey.c_str());
    internal->lodCache = NULL;
    internal->nameCache = NULL;
    internal->accessMutex = new (std::nothrow) std::mutex;
    if (!internal->accessMutex) {
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

static struct BRLObolMeshLodContext *
mesh_lod_context_create_for_db(struct db_i *dbip)
{
    if (!dbip)
	return NULL;

    char ctxName[MAXPATHLEN] = {0};
    if (dbip->dbi_filename) {
	bu_strlcpy(ctxName, dbip->dbi_filename, sizeof(ctxName));
    } else {
	snprintf(ctxName, sizeof(ctxName), "brlobol_inmem_mesh_lod_%p",
		 (void *)dbip);
    }

    return mesh_lod_context_create(ctxName);
}

static unsigned long long
mesh_lod_key_get(struct BRLObolMeshLodContext *context, const char *name)
{
    if (!context || !name)
	return 0;

    struct bu_vls keystr = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&keystr, "%s", name);
    if (bu_vls_strlen(&keystr) < 10)
	bu_vls_printf(&keystr, "GGGGGGGGGGGGG");
    unsigned long long hash = bu_data_hash(bu_vls_cstr(&keystr),
					   bu_vls_strlen(&keystr) * sizeof(char));
    bu_vls_sprintf(&keystr, "%llu:namekey", hash);

    void *data = NULL;
    std::unique_lock<std::mutex> lock(*context->i->accessMutex);
    size_t dsize = bu_cache_get(&data, bu_vls_cstr(&keystr),
				context->i->nameCache, NULL);
    lock.unlock();
    bu_vls_free(&keystr);
    if (!dsize || !data)
	return 0;

    unsigned long long meshKey = *(unsigned long long *)data;
    bu_free(data, "lod key data");
    return meshKey;
}

static int
mesh_lod_key_put(struct BRLObolMeshLodContext *context,
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

    std::unique_lock<std::mutex> lock(*context->i->accessMutex);
    size_t wsize = bu_cache_write((void *)&key, sizeof(key),
				  bu_vls_cstr(&keystr),
				  context->i->nameCache, NULL);
    lock.unlock();
    bu_vls_free(&keystr);
    return (wsize > 0) ? 0 : -1;
}

class BRLObolPopState
{
public:
    BRLObolPopState(struct BRLObolMeshLodContext *ctx,
		    const point_t *vertices,
		    size_t vertexCount,
		    const vect_t *normals,
		    const int *faces,
		    size_t faceCount,
		    unsigned long long userKey,
		    fastf_t popFaceCountThresholdRatio);
    BRLObolPopState(struct BRLObolMeshLodContext *ctx,
		    unsigned long long key);
    ~BRLObolPopState();

    int getLevel(fastf_t len);
    bool setLevel(int level);
    void shrinkMemory(void);
    void levelPoint(point_t *out, const point_t *point, int level);

    std::vector<int> lodTris;
    std::vector<fastf_t> lodTriPoints;
    std::vector<fastf_t> lodTriPointsSnapped;
    std::vector<fastf_t> lodTriNormals;

    int currLevel = -1;
    bool forceUpdate = false;
    int maxPopThresholdLevel = 0;
    bool isValid = false;
    unsigned long long hash = 0;
    point_t bbmin = VINIT_ZERO;
    point_t bbmax = VINIT_ZERO;
    struct BRLObolMeshLod *lod = NULL;

    BRLObolMeshLodDetailSetupCallback fullDetailSetup = NULL;
    BRLObolMeshLodDetailClearCallback fullDetailClear = NULL;
    BRLObolMeshLodDetailFreeCallback fullDetailFree = NULL;
    void *detailData = NULL;

private:
    void triProcess(void);
    int toLevel(int value, int level);
    fastf_t snap(fastf_t value, fastf_t min, fastf_t max, int level);
    bool isEqual(BRLObolPopRec r1, BRLObolPopRec r2, int level);
    bool triDegenerate(BRLObolPopRec r0, BRLObolPopRec r1,
		       BRLObolPopRec r2, int level);
    void cache(void);
    bool cacheTri(void);
    bool cacheWrite(const char *component, std::stringstream &stream);
    size_t cacheGet(void **data, const char *component);
    void cacheDone(void);
    bool triPopLoad(int startLevel, int level);
    void triPopTrim(int level);
    bool setupFullDetail(void);
    void clearFullDetail(void);

    float minx = FLT_MAX;
    float miny = FLT_MAX;
    float minz = FLT_MAX;
    float maxx = -FLT_MAX;
    float maxy = -FLT_MAX;
    float maxz = -FLT_MAX;
    fastf_t maxFaceRatio = 0.66;
    std::vector<unsigned short> precomputedMasks;
    size_t levelVertexCount[POP_MAXLEVEL + 1] = {0};
    size_t levelTriangleCount[POP_MAXLEVEL + 1] = {0};

    std::vector<size_t> triIndexMap;
    std::vector<size_t> vertexTriMinLevel;
    std::map<size_t, std::unordered_set<size_t>> levelTriVerts;
    std::vector<std::vector<size_t>> levelTris;

    size_t vertexCount = 0;
    const point_t *vertexArray = NULL;
    const vect_t *normalArray = NULL;
    size_t faceCount = 0;
    const int *faceArray = NULL;

    struct BRLObolMeshLodContext *context = NULL;
    struct bu_cache_txn *readTxn = NULL;
    struct bu_cache_txn *writeTxn = NULL;
    std::unique_lock<std::mutex> readLock;
};

void
BRLObolPopState::triProcess(void)
{
    vertexTriMinLevel.reserve(vertexCount);
    for (size_t vertexIndex = 0; vertexIndex < vertexCount; vertexIndex++)
	vertexTriMinLevel.push_back(POP_MAXLEVEL - 1);

    levelTris.reserve(POP_MAXLEVEL);
    for (size_t levelIndex = 0; levelIndex < POP_MAXLEVEL; levelIndex++)
	levelTris.push_back(std::vector<size_t>(0));

    for (size_t faceIndex = 0; faceIndex < faceCount; faceIndex++) {
	BRLObolPopRec triangle[3];
	bool badFace = false;
	for (size_t cornerIndex = 0; cornerIndex < 3; cornerIndex++) {
	    int faceVertex = faceArray[3 * faceIndex + cornerIndex];
	    if (faceVertex < 0 || static_cast<size_t>(faceVertex) >= vertexCount) {
		bu_log("bad face %zd - skipping\n", faceIndex);
		badFace = true;
		break;
	    }
	    triangle[cornerIndex].x = static_cast<unsigned short>(floor(
					  (vertexArray[faceVertex][X] - minx) / (maxx - minx) * USHRT_MAX));
	    triangle[cornerIndex].y = static_cast<unsigned short>(floor(
					  (vertexArray[faceVertex][Y] - miny) / (maxy - miny) * USHRT_MAX));
	    triangle[cornerIndex].z = static_cast<unsigned short>(floor(
					  (vertexArray[faceVertex][Z] - minz) / (maxz - minz) * USHRT_MAX));
	}
	if (badFace)
	    continue;

	size_t level = POP_MAXLEVEL - 1;
	for (int levelIndex = 0; levelIndex < POP_MAXLEVEL; levelIndex++) {
	    if (!triDegenerate(triangle[0], triangle[1], triangle[2],
			       levelIndex)) {
		level = static_cast<size_t>(levelIndex);
		break;
	    }
	}
	levelTris[level].push_back(faceIndex);

	for (size_t cornerIndex = 0; cornerIndex < 3; cornerIndex++) {
	    int faceVertex = faceArray[3 * faceIndex + cornerIndex];
	    if (vertexTriMinLevel[faceVertex] > level)
		vertexTriMinLevel[faceVertex] = level;
	}
    }

    for (size_t vertexIndex = 0; vertexIndex < vertexTriMinLevel.size();
	 vertexIndex++)
	levelTriVerts[vertexTriMinLevel[vertexIndex]].insert(vertexIndex);

    triIndexMap.reserve(vertexCount);
    for (size_t vertexIndex = 0; vertexIndex < vertexCount; vertexIndex++)
	triIndexMap.push_back(vertexIndex);

    size_t outputVertexIndex = 0;
    for (auto levelIt = levelTriVerts.begin(); levelIt != levelTriVerts.end();
	 ++levelIt) {
	for (auto setIt = levelIt->second.begin();
	     setIt != levelIt->second.end(); ++setIt) {
	    triIndexMap[*setIt] = outputVertexIndex;
	    outputVertexIndex++;
	}
    }

    int firstPopulatedLevel = -1;
    for (size_t levelIndex = 0; levelIndex < levelTris.size(); levelIndex++) {
	if (levelTris[levelIndex].size()) {
	    firstPopulatedLevel = static_cast<int>(levelIndex);
	    break;
	}
    }

    size_t triSum = 0;
    if (maxFaceRatio > 0.99 || maxFaceRatio < 0) {
	maxPopThresholdLevel = static_cast<int>(levelTris.size() - 1);
    } else {
	size_t faceThreshold = static_cast<size_t>(
				   static_cast<fastf_t>(faceCount) * maxFaceRatio);
	for (size_t levelIndex = 0; levelIndex < levelTris.size();
	     levelIndex++) {
	    triSum += levelTris[levelIndex].size();
	    if (triSum > faceThreshold) {
		if (triSum < faceCount)
		    maxPopThresholdLevel = static_cast<int>(levelIndex);
		else
		    maxPopThresholdLevel = levelIndex ?
					   static_cast<int>(levelIndex - 1) : 0;
		break;
	    }
	}
    }

    if (firstPopulatedLevel >= 0 && maxPopThresholdLevel < firstPopulatedLevel)
	maxPopThresholdLevel = firstPopulatedLevel;
}

BRLObolPopState::BRLObolPopState(
    struct BRLObolMeshLodContext *ctx,
    const point_t *vertices,
    size_t inputVertexCount,
    const vect_t *normals,
    const int *faces,
    size_t inputFaceCount,
    unsigned long long userKey,
    fastf_t popFaceCountThresholdRatio)
{
    context = ctx;
    maxFaceRatio = popFaceCountThresholdRatio;

    if (!userKey) {
	struct bu_data_hash_state *state = bu_data_hash_create();
	bu_data_hash_update(state, vertices, inputVertexCount * sizeof(point_t));
	bu_data_hash_update(state, faces, 3 * inputFaceCount * sizeof(int));
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

    minx = minx - fabs(MBUMP * minx);
    miny = miny - fabs(MBUMP * miny);
    minz = minz - fabs(MBUMP * minz);
    maxx = maxx + fabs(MBUMP * maxx);
    maxy = maxy + fabs(MBUMP * maxy);
    maxz = maxz + fabs(MBUMP * maxz);

    if (NEAR_EQUAL(maxx, minx, SMALL_FASTF))
	maxx = minx + 1.0f;
    if (NEAR_EQUAL(maxy, miny, SMALL_FASTF))
	maxy = miny + 1.0f;
    if (NEAR_EQUAL(maxz, minz, SMALL_FASTF))
	maxz = minz + 1.0f;

    triProcess();

    isValid = true;
    cache();
}

BRLObolPopState::BRLObolPopState(struct BRLObolMeshLodContext *ctx,
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

    {
	const char *buffer = NULL;
	size_t bufferSize = cacheGet((void **)&buffer, CACHE_POP_MAX_LEVEL);
	if (!bufferSize || bufferSize != sizeof(maxPopThresholdLevel)) {
	    cacheDone();
	    return;
	}
	memcpy(&maxPopThresholdLevel, buffer, sizeof(maxPopThresholdLevel));
	cacheDone();
    }

    {
	const char *buffer = NULL;
	size_t bufferSize = cacheGet((void **)&buffer, CACHE_POP_SWITCH_LEVEL);
	if (bufferSize && bufferSize != sizeof(maxFaceRatio)) {
	    cacheDone();
	    return;
	}
	if (bufferSize)
	    memcpy(&maxFaceRatio, buffer, sizeof(maxFaceRatio));
	else
	    maxFaceRatio = 0.66;
	cacheDone();
    }

    {
	const char *buffer = NULL;
	size_t bufferSize = cacheGet((void **)&buffer, CACHE_VERTEX_COUNT);
	if (bufferSize != sizeof(levelVertexCount)) {
	    cacheDone();
	    return;
	}
	memcpy(&levelVertexCount, buffer, sizeof(levelVertexCount));
	cacheDone();
    }

    {
	const char *buffer = NULL;
	size_t bufferSize = cacheGet((void **)&buffer, CACHE_TRI_COUNT);
	if (bufferSize != sizeof(levelTriangleCount)) {
	    cacheDone();
	    return;
	}
	memcpy(&levelTriangleCount, buffer, sizeof(levelTriangleCount));
	cacheDone();
    }

    {
	float minmax[6];
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
	cacheDone();
    }

    if (!setLevel(0))
	return;
    isValid = true;
}

BRLObolPopState::~BRLObolPopState()
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

bool
BRLObolPopState::triPopLoad(int startLevel, int level)
{
    struct bu_vls keyBuffer = BU_VLS_INIT_ZERO;

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
	cacheDone();
    }

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
	cacheDone();
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
	cacheDone();
    }

    bu_vls_free(&keyBuffer);
    return true;
}

void
BRLObolPopState::shrinkMemory(void)
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
BRLObolPopState::triPopTrim(int level)
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

int
BRLObolPopState::getLevel(fastf_t viewLength)
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

bool
BRLObolPopState::setupFullDetail(void)
{
    if (!lod || !fullDetailSetup)
	return false;

    struct BRLObolMeshLodDetail detail;
    brlobol_mesh_lod_detail_init(&detail);
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
BRLObolPopState::clearFullDetail(void)
{
    if (fullDetailClear)
	fullDetailClear(detailData);
    mesh_lod_active_data_clear(lod);
}

bool
BRLObolPopState::setLevel(int level)
{
    if (level > maxPopThresholdLevel && !fullDetailSetup)
	level = maxPopThresholdLevel;

    if (level == currLevel && !forceUpdate)
	return true;

    if (forceUpdate) {
	forceUpdate = false;
	clearFullDetail();
	shrinkMemory();
	currLevel = -1;
	return setLevel(level);
    }

    if (level > currLevel && level <= maxPopThresholdLevel) {
	if (!lodTriPoints.size()) {
	    if (!triPopLoad(-1, level))
		return false;
	} else if (!triPopLoad(currLevel, level)) {
	    return false;
	}
    }

    if (level < currLevel && level <= maxPopThresholdLevel &&
	currLevel <= maxPopThresholdLevel) {
	if (!lodTriPoints.size()) {
	    if (!triPopLoad(-1, level))
		return false;
	} else {
	    triPopTrim(level);
	}
    }

    if (level < currLevel && level <= maxPopThresholdLevel &&
	currLevel > maxPopThresholdLevel) {
	clearFullDetail();
	if (!triPopLoad(-1, level))
	    return false;
    }

    if (level > currLevel && level > maxPopThresholdLevel &&
	currLevel <= maxPopThresholdLevel) {
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
BRLObolPopState::cacheWrite(const char *component, std::stringstream &stream)
{
    std::string buffer = stream.str();
    if (!buffer.length())
	return false;

    std::string keystr = std::to_string(hash) + std::string(":") +
			 std::string(component);
    size_t wsize = bu_cache_write((void *)buffer.data(), buffer.length(),
				  keystr.c_str(), context->i->lodCache,
				  &writeTxn);
    return (wsize > 0);
}

size_t
BRLObolPopState::cacheGet(void **data, const char *component)
{
    if (context && context->i && context->i->accessMutex &&
	!readLock.owns_lock())
	readLock = std::unique_lock<std::mutex>(*context->i->accessMutex);
    std::string keystr = std::to_string(hash) + std::string(":") +
			 std::string(component);
    return bu_cache_get(data, keystr.c_str(), context->i->lodCache, &readTxn);
}

void
BRLObolPopState::cacheDone(void)
{
    bu_cache_get_done(&readTxn);
    if (readLock.owns_lock())
	readLock.unlock();
}

bool
BRLObolPopState::cacheTri(void)
{
    {
	std::stringstream stream;
	stream.write(reinterpret_cast<const char *>(&maxPopThresholdLevel),
		     sizeof(maxPopThresholdLevel));
	if (!cacheWrite(CACHE_POP_MAX_LEVEL, stream))
	    return false;
    }

    {
	std::stringstream stream;
	stream.write(reinterpret_cast<const char *>(&maxFaceRatio),
		     sizeof(maxFaceRatio));
	if (!cacheWrite(CACHE_POP_SWITCH_LEVEL, stream))
	    return false;
    }

    {
	std::stringstream stream;
	for (size_t levelIndex = 0; levelIndex <= POP_MAXLEVEL; levelIndex++) {
	    size_t count = 0;
	    if (levelTriVerts.find(levelIndex) == levelTriVerts.end()) {
		stream.write(reinterpret_cast<const char *>(&count),
			     sizeof(count));
		continue;
	    }
	    if (static_cast<int>(levelIndex) > maxPopThresholdLevel ||
		!levelTriVerts[levelIndex].size()) {
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
	    if (static_cast<int>(levelIndex) > maxPopThresholdLevel ||
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

    for (int levelIndex = 0; levelIndex <= maxPopThresholdLevel;
	 levelIndex++) {
	std::stringstream stream;
	if (levelTriVerts.find(levelIndex) == levelTriVerts.end())
	    continue;
	if (!levelTriVerts[levelIndex].size())
	    continue;
	for (auto setIt = levelTriVerts[levelIndex].begin();
	     setIt != levelTriVerts[levelIndex].end(); ++setIt) {
	    point_t point;
	    VMOVE(point, vertexArray[*setIt]);
	    stream.write(reinterpret_cast<const char *>(&point[0]), sizeof(point_t));
	}
	bu_vls_sprintf(&keyBuffer, "%s%d", CACHE_VERT_LEVEL, levelIndex);
	if (!cacheWrite(bu_vls_cstr(&keyBuffer), stream)) {
	    bu_vls_free(&keyBuffer);
	    return false;
	}
    }

    for (int levelIndex = 0; levelIndex <= maxPopThresholdLevel;
	 levelIndex++) {
	std::stringstream stream;
	if (!levelTris[levelIndex].size())
	    continue;
	for (auto triIt = levelTris[levelIndex].begin();
	     triIt != levelTris[levelIndex].end(); ++triIt) {
	    int vertices[3];
	    vertices[0] = static_cast<int>(triIndexMap[faceArray[3 * (*triIt) + 0]]);
	    vertices[1] = static_cast<int>(triIndexMap[faceArray[3 * (*triIt) + 1]]);
	    vertices[2] = static_cast<int>(triIndexMap[faceArray[3 * (*triIt) + 2]]);
	    stream.write(reinterpret_cast<const char *>(&vertices[0]),
			 sizeof(vertices));
	}
	bu_vls_sprintf(&keyBuffer, "%s%d", CACHE_TRI_LEVEL, levelIndex);
	if (!cacheWrite(bu_vls_cstr(&keyBuffer), stream)) {
	    bu_vls_free(&keyBuffer);
	    return false;
	}
    }

    if (normalArray) {
	for (int levelIndex = 0; levelIndex <= maxPopThresholdLevel;
	     levelIndex++) {
	    std::stringstream stream;
	    if (!levelTris[levelIndex].size())
		continue;
	    for (auto triIt = levelTris[levelIndex].begin();
		 triIt != levelTris[levelIndex].end(); ++triIt) {
		for (int cornerIndex = 0; cornerIndex < 3; cornerIndex++) {
		    vect_t normal;
		    const int normalIndex = 3 * (*triIt) + cornerIndex;
		    VMOVE(normal, normalArray[normalIndex]);
		    stream.write(reinterpret_cast<const char *>(&normal[0]),
				 sizeof(vect_t));
		}
	    }
	    bu_vls_sprintf(&keyBuffer, "%s%d", CACHE_VERTNORM_LEVEL,
			   levelIndex);
	    if (!cacheWrite(bu_vls_cstr(&keyBuffer), stream)) {
		bu_vls_free(&keyBuffer);
		return false;
	    }
	}
    }

    bu_vls_free(&keyBuffer);
    return true;
}

void
BRLObolPopState::cache(void)
{
    if (!hash) {
	isValid = false;
	return;
    }

    int writeSem = mesh_lod_cache_write_semaphore();
    bu_semaphore_acquire(writeSem);
    std::unique_lock<std::mutex> cacheLock(*context->i->accessMutex);

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
BRLObolPopState::toLevel(int value, int level)
{
    return static_cast<int>(floor(value / double(precomputedMasks[level])));
}

fastf_t
BRLObolPopState::snap(fastf_t value, fastf_t min, fastf_t max, int level)
{
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
BRLObolPopState::levelPoint(point_t *out, const point_t *point, int level)
{
    fastf_t nx = snap((*point)[X], minx, maxx, level);
    fastf_t ny = snap((*point)[Y], miny, maxy, level);
    fastf_t nz = snap((*point)[Z], minz, maxz, level);
    VSET(*out, nx, ny, nz);
}

bool
BRLObolPopState::isEqual(BRLObolPopRec r1, BRLObolPopRec r2, int level)
{
    bool equalX = (toLevel(r1.x, level) == toLevel(r2.x, level));
    bool equalY = (toLevel(r1.y, level) == toLevel(r2.y, level));
    bool equalZ = (toLevel(r1.z, level) == toLevel(r2.z, level));
    return (equalX && equalY && equalZ);
}

bool
BRLObolPopState::triDegenerate(BRLObolPopRec r0,
			       BRLObolPopRec r1,
			       BRLObolPopRec r2,
			       int level)
{
    return isEqual(r0, r1, level) || isEqual(r1, r2, level) ||
	   isEqual(r0, r2, level);
}

static void
mesh_lod_active_pop_data_publish(struct BRLObolMeshLod *lod,
				 BRLObolPopState *state)
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
    lod->points = reinterpret_cast<const point_t *>(
		      state->lodTriPointsSnapped.data());
    lod->pcnt = static_cast<int>(state->lodTriPointsSnapped.size() / 3);
}

static unsigned long long
mesh_lod_cache_generate(struct BRLObolMeshLodContext *context,
			const point_t *vertices,
			size_t vertexCount,
			const vect_t *normals,
			const int *faces,
			size_t faceCount,
			unsigned long long userKey,
			fastf_t fidelityRatio)
{
    if (!context || !vertices || !vertexCount || !faces || !faceCount)
	return 0;

    BRLObolPopState state(context, vertices, vertexCount, normals, faces,
			  faceCount, userKey, fidelityRatio);
    return state.isValid ? state.hash : 0;
}

static struct BRLObolMeshLod *
mesh_lod_create(struct BRLObolMeshLodContext *context,
		unsigned long long key)
{
    if (!context || !key)
	return NULL;

    BRLObolPopState *state = new (std::nothrow) BRLObolPopState(context, key);
    if (!state)
	return NULL;

    if (!state->isValid) {
	delete state;
	return NULL;
    }

    struct BRLObolMeshLod *lod;
    BU_GET(lod, struct BRLObolMeshLod);
    lod->context = context;
    lod->state = state;
    mesh_lod_active_data_clear(lod);
    VMOVE(lod->bmin, state->bbmin);
    VMOVE(lod->bmax, state->bbmax);
    state->lod = lod;
    return lod;
}

static int
mesh_lod_level(struct BRLObolMeshLod *lod, int level, int reset)
{
    if (!lod || !lod->state)
	return -1;
    BRLObolPopState *state = lod->state;
    if (level < 0)
	return state->currLevel;

    if (reset)
	state->forceUpdate = true;
    if (!state->setLevel(level)) {
	mesh_lod_active_data_clear(lod);
	return -1;
    }

    if (state->currLevel <= state->maxPopThresholdLevel)
	mesh_lod_active_pop_data_publish(lod, state);

    return state->currLevel;
}

static int
mesh_lod_payload_available(struct BRLObolMeshLodContext *context,
			   unsigned long long key)
{
    if (!context || !key)
	return 0;

    struct BRLObolMeshLod *lod = mesh_lod_create(context, key);
    if (!lod)
	return 0;

    lod->context = NULL;
    brlobol_mesh_lod_destroy(lod);
    return 1;
}

static void
mesh_lod_status_current(struct db_i *dbip,
			struct BRLObolMeshLodContext *context,
			const char *name,
			struct BRLObolMeshLodCacheStatus *status)
{
    if (!status)
	return;

    brlobol_mesh_lod_cache_status_init(status);
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
mesh_lod_cache_del(struct BRLObolMeshLodContext *context,
		   unsigned long long hash,
		   const char *component)
{
    std::string keystr = std::to_string(hash) + std::string(":") +
			 std::string(component);
    bu_cache_clear(keystr.c_str(), context->i->lodCache, NULL);
}

static void
mesh_lod_cache_clear_context(struct BRLObolMeshLodContext *context,
			     unsigned long long key)
{
    char dir[MAXPATHLEN];
    std::unique_lock<std::mutex> lock;

    if (context && context->i && context->i->accessMutex)
	lock = std::unique_lock<std::mutex>(*context->i->accessMutex);

    if (context && key) {
	mesh_lod_cache_del(context, key, CACHE_POP_MAX_LEVEL);
	mesh_lod_cache_del(context, key, CACHE_POP_SWITCH_LEVEL);
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
	return;
    }

    if (context && !key) {
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
mesh_lod_bot_detail_setup(struct BRLObolMeshLodDetail *detail, void *cbData)
{
    if (!detail || !cbData)
	return -1;

    struct BRLObolMeshLodBotDetail *callbackData =
	    static_cast<struct BRLObolMeshLodBotDetail *>(cbData);
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
    struct BRLObolMeshLodBotDetail *callbackData =
	    static_cast<struct BRLObolMeshLodBotDetail *>(cbData);

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
	struct BRLObolMeshLodBotDetail *callbackData =
		static_cast<struct BRLObolMeshLodBotDetail *>(cbData);
	BU_PUT(callbackData, struct BRLObolMeshLodBotDetail);
    }

    return 0;
}

void
brlobol_mesh_lod_cache_init(struct db_i *dbip, int verbose)
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

    struct BRLObolMeshLodCacheStatus status =
	    BRLOBOL_MESH_LOD_CACHE_STATUS_INIT;
    if (brlobol_mesh_lod_cache_status(dbip, dp->d_namep, &status) == BRLCAD_OK &&
	status.has_cache_key && status.has_cached_payload)
	continue;

    if (verbose > 1)
	bu_log("Processing(%d):  %s\n", completed + 1, dp->d_namep);

    (void)brlobol_mesh_lod_cache_update(dbip, dp->d_namep);
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
brlobol_mesh_lod_cache_clear_database(struct db_i *dbip)
{
    struct BRLObolMeshLodContext *context = mesh_lod_context_create_for_db(dbip);
    if (!context)
	return;
    mesh_lod_cache_clear_context(context, 0);
    mesh_lod_context_destroy(context);
}

void
brlobol_mesh_lod_cache_clear_all(void)
{
    mesh_lod_cache_clear_context(NULL, 0);
}

int
brlobol_mesh_lod_cache_update(struct db_i *dbip, const char *name)
{
    return brlobol_mesh_lod_cache_refresh(dbip, name, NULL);
}

void
brlobol_mesh_lod_cache_status_init(struct BRLObolMeshLodCacheStatus *status)
{
    struct BRLObolMeshLodCacheStatus defaults =
	    BRLOBOL_MESH_LOD_CACHE_STATUS_INIT;
    if (status)
	*status = defaults;
}

int
brlobol_mesh_lod_cache_status(struct db_i *dbip,
			      const char *name,
			      struct BRLObolMeshLodCacheStatus *status)
{
    if (!dbip || !name || !status)
	return BRLCAD_ERROR;

    struct BRLObolMeshLodContext *context = mesh_lod_context_create_for_db(dbip);
    if (!context)
	return BRLCAD_ERROR;

    mesh_lod_status_current(dbip, context, name, status);
    mesh_lod_context_destroy(context);
    return BRLCAD_OK;
}

int
brlobol_mesh_lod_cache_invalidate(struct db_i *dbip,
				  const char *name,
				  struct BRLObolMeshLodCacheStatus *status)
{
    struct BRLObolMeshLodCacheStatus current =
	    BRLOBOL_MESH_LOD_CACHE_STATUS_INIT;

    if (status)
	brlobol_mesh_lod_cache_status_init(status);
    if (!dbip || !name)
	return BRLCAD_ERROR;

    struct BRLObolMeshLodContext *context = mesh_lod_context_create_for_db(dbip);
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
brlobol_mesh_lod_cache_refresh(struct db_i *dbip,
			       const char *name,
			       struct BRLObolMeshLodCacheStatus *status)
{
    struct BRLObolMeshLodCacheStatus current =
	    BRLOBOL_MESH_LOD_CACHE_STATUS_INIT;

    if (status)
	brlobol_mesh_lod_cache_status_init(status);
    if (!dbip || !name)
	return BRLCAD_ERROR;

    struct BRLObolMeshLodContext *context = mesh_lod_context_create_for_db(dbip);
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

    std::vector<fastf_t> generatedNormals;
    const point_t *botVertices =
	reinterpret_cast<const point_t *>(bot->vertices);
    mesh_lod_generate_crease_normals(generatedNormals, botVertices,
	bot->num_vertices, bot->faces, bot->num_faces);
    const vect_t *botNormals = generatedNormals.empty() ? NULL :
	reinterpret_cast<const vect_t *>(generatedNormals.data());

    unsigned long long key = mesh_lod_cache_generate(
				 context, botVertices, bot->num_vertices,
				 botNormals, bot->faces, bot->num_faces, 0, 0.66);
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
brlobol_mesh_lod_cache_store_mesh(
    struct db_i *dbip,
    const char *name,
    const point_t *vertices,
    size_t vertexCount,
    const vect_t *normals,
    const int *faces,
    size_t faceCount,
    unsigned long long userKey,
    fastf_t fidelityRatio,
    struct BRLObolMeshLodCacheStatus *status)
{
    struct BRLObolMeshLodCacheStatus current =
	    BRLOBOL_MESH_LOD_CACHE_STATUS_INIT;

    if (status)
	brlobol_mesh_lod_cache_status_init(status);
    if (!dbip || !name)
	return BRLCAD_ERROR;

    struct BRLObolMeshLodContext *context = mesh_lod_context_create_for_db(dbip);
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

    std::vector<fastf_t> generatedNormals;
    if (!normals)
	mesh_lod_generate_crease_normals(generatedNormals, vertices,
	    vertexCount, faces, faceCount);
    const vect_t *cacheNormals = normals ? normals :
	(generatedNormals.empty() ? NULL :
	 reinterpret_cast<const vect_t *>(generatedNormals.data()));

    unsigned long long key = mesh_lod_cache_generate(
				 context, vertices, vertexCount, cacheNormals, faces, faceCount, userKey,
				 fidelityRatio);
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

struct BRLObolMeshLod *
brlobol_mesh_lod_get(struct db_i *dbip, const char *name)
{
    if (!dbip || !name)
	return NULL;

    struct BRLObolMeshLodContext *context = mesh_lod_context_create_for_db(dbip);
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

    struct BRLObolMeshLod *lod = mesh_lod_create(context, key);
    if (!lod) {
	mesh_lod_context_destroy(context);
	return NULL;
    }

    if (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) {
	struct BRLObolMeshLodBotDetail *callbackData;
	BU_GET(callbackData, struct BRLObolMeshLodBotDetail);
	memset(callbackData, 0, sizeof(*callbackData));
	callbackData->dbip = dbip;
	callbackData->dp = dp;
	if (!brlobol_mesh_lod_detail_callbacks_set(
		lod, mesh_lod_bot_detail_setup, mesh_lod_bot_detail_clear,
		mesh_lod_bot_detail_free, callbackData))
	    mesh_lod_bot_detail_free(callbackData);
    }

    return lod;
}

int
brlobol_mesh_lod_load_level(struct BRLObolMeshLod *lod, int level, int reset)
{
    return mesh_lod_level(lod, level, reset);
}

int
brlobol_mesh_lod_load_view(struct BRLObolMeshLod *lod,
			   const struct rt_view_info *info,
			   int reset)
{
    if (!lod || !lod->state)
	return -1;

    struct rt_view_info sanitized = RT_VIEW_INFO_INIT;
    if (info)
	sanitized = *info;
    rt_view_info_sanitize(&sanitized);

    fastf_t policyScale = sanitized.lod.scale;
    if (sanitized.size <= SMALL_FASTF)
	sanitized.size = 1.0;
    if (policyScale <= SMALL_FASTF)
	policyScale = 1.0;

    BRLObolPopState *state = lod->state;
    int viewLevel = static_cast<int>(
			static_cast<double>(state->getLevel(sanitized.size)) * policyScale);
    viewLevel = (viewLevel < 0) ? 0 : viewLevel;
    viewLevel = (viewLevel >= POP_MAXLEVEL) ? POP_MAXLEVEL - 1 : viewLevel;

    return mesh_lod_level(lod, viewLevel, reset);
}

int
brlobol_mesh_lod_load_view_scene_ref(struct BRLObolMeshLod *lod,
				     rt_view_scene_ref visibility_ref,
				     void *view_ctx,
				     int reset)
{
    (void)visibility_ref;
    struct rt_view_info info = RT_VIEW_INFO_INIT;
    if (view_ctx)
	rt_view_context_info_get(&info, view_ctx);
    return brlobol_mesh_lod_load_view(lod, &info, reset);
}

void
brlobol_mesh_lod_free_scene_ref(rt_view_scene_ref ref)
{
    (void)ref;
}

int
brlobol_mesh_lod_current_level(const struct BRLObolMeshLod *lod)
{
    return (lod && lod->state) ? lod->state->currLevel : -1;
}

int
brlobol_mesh_lod_has_active_data(const struct BRLObolMeshLod *lod)
{
    return (lod && lod->faces && lod->fcnt > 0 &&
	    lod->points && lod->pcnt > 0 &&
	    lod->points_orig && lod->porig_cnt > 0) ? 1 : 0;
}

int
brlobol_mesh_lod_data_get(const struct BRLObolMeshLod *lod,
			  struct BRLObolMeshLodData *data)
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

    return brlobol_mesh_lod_has_active_data(lod);
}

void
brlobol_mesh_lod_info_init(struct BRLObolMeshLodInfo *info)
{
    struct BRLObolMeshLodInfo defaults = BRLOBOL_MESH_LOD_INFO_INIT;
    if (info)
	*info = defaults;
}

int
brlobol_mesh_lod_info_get(const struct BRLObolMeshLod *lod,
			  struct BRLObolMeshLodInfo *info)
{
    if (!info)
	return 0;

    brlobol_mesh_lod_info_init(info);
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
    VMOVE(info->bmin, lod->bmin);
    VMOVE(info->bmax, lod->bmax);

    return info->has_faces && info->has_points && info->has_original_points;
}

void
brlobol_mesh_lod_detail_init(struct BRLObolMeshLodDetail *detail)
{
    if (detail)
	memset(detail, 0, sizeof(*detail));
}

int
brlobol_mesh_lod_detail_callbacks_set(
    struct BRLObolMeshLod *lod,
    BRLObolMeshLodDetailSetupCallback setupCallback,
    BRLObolMeshLodDetailClearCallback clearCallback,
    BRLObolMeshLodDetailFreeCallback freeCallback,
    void *cbData)
{
    if (!lod || !lod->state || !setupCallback)
	return 0;

    BRLObolPopState *state = lod->state;
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
	} else if (state->currLevel > state->maxPopThresholdLevel) {
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
brlobol_mesh_lod_detail_callbacks_clear(struct BRLObolMeshLod *lod)
{
    if (!lod || !lod->state)
	return;

    BRLObolPopState *state = lod->state;
    if (state->fullDetailFree)
	state->fullDetailFree(state->detailData);
    else if (state->fullDetailClear)
	state->fullDetailClear(state->detailData);
    state->fullDetailSetup = NULL;
    state->fullDetailClear = NULL;
    state->fullDetailFree = NULL;
    state->detailData = NULL;
    if (state->currLevel > state->maxPopThresholdLevel) {
	state->currLevel = -1;
	mesh_lod_active_data_clear(lod);
    }
}

void
brlobol_mesh_lod_memshrink(struct BRLObolMeshLod *lod)
{
    if (!lod || !lod->state)
	return;

    if (lod->state->currLevel > lod->state->maxPopThresholdLevel)
	return;

    lod->state->shrinkMemory();
    lod->state->forceUpdate = true;
    mesh_lod_active_data_clear(lod);
}

void
brlobol_mesh_lod_destroy(struct BRLObolMeshLod *lod)
{
    if (!lod)
	return;

    BRLObolMeshLodContext *context = lod->context;
    delete lod->state;
    lod->state = NULL;
    lod->context = NULL;
    BU_PUT(lod, struct BRLObolMeshLod);
    mesh_lod_context_destroy(context);
}
