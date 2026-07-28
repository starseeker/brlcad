/*                 G E N E R A T E _ U N I Q U E _ M E S H _ S T R E S S . C P P
 * BRL-CAD
 *
 * Generate a vehicle-like workload with many genuinely distinct, moderately
 * detailed BoT leaves.  This complements repeated-instance Lucy stress:
 * every leaf has different vertex content and therefore exercises distinct
 * import, PoP construction, cache, publication, upload, and residency work.
 */

#include "common.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

#include "bn/mat.h"
#include "raytrace.h"
#include "wdb.h"

namespace {

bool
parse_positive(const char *text, size_t maximum, size_t &value)
{
    if (!text || !text[0])
	return false;
    char *end = NULL;
    const unsigned long long parsed = std::strtoull(text, &end, 10);
    if (!end || *end || parsed == 0 || parsed > maximum)
	return false;
    value = static_cast<size_t>(parsed);
    return true;
}

size_t
mesh_grid_cells(size_t mesh, size_t maximum)
{
    /* A deliberately skewed distribution spans approximately 500 through
     * 100,000 faces while keeping a 5,000-leaf aggregate practical.  The
     * source-size bucket is independent of physical extent and placement. */
    const size_t percentile = mesh % 100;
    size_t cells = 224;          /* 100,352 faces */
    if (percentile < 25)
	cells = 16;               /* 512 faces */
    else if (percentile < 45)
	cells = 24;               /* 1,152 faces */
    else if (percentile < 60)
	cells = 32;               /* 2,048 faces */
    else if (percentile < 72)
	cells = 48;               /* 4,608 faces */
    else if (percentile < 82)
	cells = 64;               /* 8,192 faces */
    else if (percentile < 90)
	cells = 96;               /* 18,432 faces */
    else if (percentile < 95)
	cells = 128;              /* 32,768 faces */
    else if (percentile < 98)
	cells = 160;              /* 51,200 faces */
    else if (percentile < 99)
	cells = 192;              /* 73,728 faces */
    return std::max<size_t>(1, std::min(cells, maximum));
}

bool
finish_group(struct rt_wdb *wdbp, struct wmember &members,
    std::vector<std::string> &groupNames, size_t groupIndex)
{
    if (BU_LIST_IS_EMPTY(&members.l))
	return true;
    const bool region = groupIndex % 4 != 3;
    char name[64] = {0};
    std::snprintf(name, sizeof(name), region ?
	"unique_region_%06zu.r" : "unique_leaf_group_%06zu.c", groupIndex);
    unsigned char color[3] = {
	static_cast<unsigned char>(48 + (groupIndex * 53) % 192),
	static_cast<unsigned char>(48 + (groupIndex * 97) % 192),
	static_cast<unsigned char>(48 + (groupIndex * 151) % 192)
    };
    if (mk_lcomb(wdbp, name, &members, region ? 1 : 0,
	    NULL, NULL, color, 0) != 0)
	return false;
    groupNames.emplace_back(name);
    BU_LIST_INIT(&members.l);
    return true;
}

bool
finish_hierarchy_level(struct rt_wdb *wdbp,
    const std::vector<std::string> &children,
    std::vector<std::string> &parents, size_t level)
{
    static const size_t childrenPerParent = 8;
    for (size_t first = 0; first < children.size();
	 first += childrenPerParent) {
	struct wmember members;
	BU_LIST_INIT(&members.l);
	const size_t last =
	    std::min(children.size(), first + childrenPerParent);
	for (size_t child = first; child < last; ++child) {
	    if (!mk_addmember(
		    children[child].c_str(), &members.l, NULL, WMOP_UNION))
		return false;
	}
	char name[64] = {0};
	std::snprintf(name, sizeof(name), "unique_level_%02zu_%06zu.c",
	    level, first / childrenPerParent);
	if (mk_lcomb(wdbp, name, &members, 0,
		NULL, NULL, NULL, 0) != 0)
	    return false;
	parents.emplace_back(name);
    }
    return true;
}

} // namespace

int
main(int argc, char **argv)
{
    if (argc < 2 || argc > 4) {
	std::fprintf(stderr,
	    "Usage: %s output.g [mesh-count=5000] [max-grid-cells=224]\n",
	    argv[0]);
	return 2;
    }

    size_t meshCount = 5000;
    size_t maxGridCells = 224;
    if ((argc > 2 && !parse_positive(argv[2], 10000, meshCount)) ||
	(argc > 3 && !parse_positive(argv[3], 256, maxGridCells))) {
	std::fprintf(stderr, "mesh count or grid cell count is out of range\n");
	return 2;
    }

    struct db_i *dbip = db_create(argv[1], 5);
    if (!dbip) {
	std::fprintf(stderr, "cannot create %s\n", argv[1]);
	return 1;
    }
    struct rt_wdb *wdbp = wdb_dbopen(dbip, RT_WDB_TYPE_DB_DISK);
    if (!wdbp) {
	db_close(dbip);
	std::fprintf(stderr, "cannot open %s for writing\n", argv[1]);
	return 1;
    }
    mk_id(wdbp, "qged distinct detailed mesh LoD stress");

    static const size_t leavesPerGroup = 16;
    struct wmember groupMembers;
    BU_LIST_INIT(&groupMembers.l);
    std::vector<std::string> groupNames;
    /* Fill an elongated three-dimensional envelope instead of laying every
     * panel in one plane.  A coplanar fixture becomes nearly edge-on at
     * ae 90 0 and accidentally tests tiny projected cuts rather than the
     * intended many-visible-detailed-mesh workload. */
    const size_t columns = static_cast<size_t>(std::ceil(std::cbrt(
	static_cast<double>(meshCount) * 4.0)));
    const size_t rows = static_cast<size_t>(std::ceil(std::sqrt(
	static_cast<double>(meshCount) /
	static_cast<double>(std::max<size_t>(1, columns)))));
    const size_t layerCapacity = std::max<size_t>(1, columns * rows);
    const size_t layers =
	(meshCount + layerCapacity - 1) / layerCapacity;
    const double spacing = 218.0;
    bool ok = true;
    size_t aggregateFaces = 0;
    size_t minimumFaces = std::numeric_limits<size_t>::max();
    size_t maximumFaces = 0;

    for (size_t mesh = 0; mesh < meshCount && ok; ++mesh) {
	const size_t gridCells = mesh_grid_cells(mesh, maxGridCells);
	const size_t side = gridCells + 1;
	const size_t vertexCount = side * side;
	const size_t faceCount = gridCells * gridCells * 2;
	std::vector<fastf_t> vertices(vertexCount * 3);
	std::vector<int> faces(faceCount * 3);
	size_t faceOffset = 0;
	for (size_t y = 0; y < gridCells; ++y) {
	    for (size_t x = 0; x < gridCells; ++x) {
		const int v00 = static_cast<int>(y * side + x);
		const int v10 = v00 + 1;
		const int v01 = static_cast<int>((y + 1) * side + x);
		const int v11 = v01 + 1;
		faces[faceOffset++] = v00;
		faces[faceOffset++] = v10;
		faces[faceOffset++] = v11;
		faces[faceOffset++] = v00;
		faces[faceOffset++] = v11;
		faces[faceOffset++] = v01;
	    }
	}
	const double phase = static_cast<double>(mesh) * 0.61803398875;
	const double amplitude = 2.0 + static_cast<double>(mesh % 17) * 0.11;
	const double meshSpan =
	    45.0 + static_cast<double>((mesh * 73) % 146);
	const size_t featureX = (mesh * 17 + 3) % side;
	const size_t featureY = (mesh * 29 + 5) % side;
	for (size_t y = 0; y < side; ++y) {
	    for (size_t x = 0; x < side; ++x) {
		const size_t vertex = y * side + x;
		const double nx =
		    static_cast<double>(x) / static_cast<double>(gridCells);
		const double ny =
		    static_cast<double>(y) / static_cast<double>(gridCells);
		const double dx = static_cast<double>(x) -
		    static_cast<double>(featureX);
		const double dy = static_cast<double>(y) -
		    static_cast<double>(featureY);
		const double bump = std::exp(
		    -(dx * dx + dy * dy) /
		    std::max(4.0, static_cast<double>(gridCells) * 0.35));
		vertices[vertex * 3 + X] = (nx - 0.5) * meshSpan;
		vertices[vertex * 3 + Y] = (ny - 0.5) * meshSpan;
		vertices[vertex * 3 + Z] =
		    amplitude * std::sin(nx * 11.0 + phase) *
			std::cos(ny * 9.0 - phase * 0.7) +
		    bump * (3.0 + static_cast<double>(mesh % 13) * 0.2);
	    }
	}

	char leafName[64] = {0};
	std::snprintf(leafName, sizeof(leafName),
	    "unique_mesh_%06zu.bot", mesh);
	if (mk_bot(wdbp, leafName, RT_BOT_SURFACE, RT_BOT_CCW, 0,
		static_cast<int>(vertexCount), static_cast<int>(faceCount),
		vertices.data(), faces.data(), NULL, NULL) != 0) {
	    std::fprintf(stderr, "failed writing %s\n", leafName);
	    ok = false;
	    break;
	}
	aggregateFaces += faceCount;
	minimumFaces = std::min(minimumFaces, faceCount);
	maximumFaces = std::max(maximumFaces, faceCount);

	mat_t placement;
	const size_t column = mesh % columns;
	const size_t row = (mesh / columns) % rows;
	const size_t layer = mesh / layerCapacity;
	const double tx =
	    (static_cast<double>(column) -
		(static_cast<double>(columns) - 1.0) * 0.5) * spacing;
	const double ty =
	    (static_cast<double>(row) -
		(static_cast<double>(rows) - 1.0) * 0.5) * spacing;
	const double tz =
	    (static_cast<double>(layer) -
		(static_cast<double>(layers) - 1.0) * 0.5) * spacing +
	    18.0 * std::sin(static_cast<double>(column) * 0.47) +
	    12.0 * std::cos(static_cast<double>(row) * 0.53);
	/* Exercise projection, PCA, lighting, and culling with panels spanning
	 * all principal orientations plus a deterministic twist. */
	const double twist = std::fmod(
	    static_cast<double>(mesh) * 137.50776405, 360.0);
	double rotateX = 7.0 * std::sin(phase);
	double rotateY = 7.0 * std::cos(phase * 0.7);
	switch (mesh % 6) {
	    case 1: rotateX += 90.0; break;
	    case 2: rotateX -= 90.0; break;
	    case 3: rotateY += 90.0; break;
	    case 4: rotateY -= 90.0; break;
	    case 5:
		rotateX += 45.0;
		rotateY += 35.0;
		break;
	    default: break;
	}
	bn_mat_angles(placement, rotateX, rotateY, twist);
	MAT_DELTAS(placement, tx, ty, tz);
	if (!mk_addmember(
		leafName, &groupMembers.l, placement, WMOP_UNION)) {
	    std::fprintf(stderr, "failed grouping %s\n", leafName);
	    ok = false;
	    break;
	}
	if ((mesh + 1) % leavesPerGroup == 0 &&
	    !finish_group(wdbp, groupMembers, groupNames,
		groupNames.size()))
	    ok = false;
    }
    if (ok && !finish_group(
	    wdbp, groupMembers, groupNames, groupNames.size()))
	ok = false;

    std::vector<std::string> hierarchy = groupNames;
    size_t hierarchyLevel = 0;
    while (ok && hierarchy.size() > 1) {
	std::vector<std::string> parents;
	if (!finish_hierarchy_level(
		wdbp, hierarchy, parents, hierarchyLevel++)) {
	    ok = false;
	    break;
	}
	hierarchy.swap(parents);
    }
    if (ok) {
	struct wmember topMembers;
	BU_LIST_INIT(&topMembers.l);
	for (const std::string &root : hierarchy) {
	    if (!mk_addmember(
		    root.c_str(), &topMembers.l, NULL, WMOP_UNION)) {
		ok = false;
		break;
	    }
	}
	if (ok && mk_lcomb(wdbp, "unique_mesh_stress", &topMembers, 0,
		NULL, NULL, NULL, 0) != 0)
	    ok = false;
    }

    wdb_close(wdbp);
    if (!ok)
	return 1;
    std::printf("generated %zu distinct meshes, %zu..%zu faces each, "
	"%zu aggregate faces in %s\n",
	meshCount, minimumFaces == std::numeric_limits<size_t>::max() ?
	    0 : minimumFaces,
	maximumFaces, aggregateFaces, argv[1]);
    return 0;
}
