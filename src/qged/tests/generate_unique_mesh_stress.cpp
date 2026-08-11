/*                 G E N E R A T E _ U N I Q U E _ M E S H _ S T R E S S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
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
#include "bu/app.h"
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

double
mesh_physical_unit(size_t mesh)
{
    const uint64_t hash = static_cast<uint64_t>(mesh) *
	11400714819323198485ULL + 0x9e3779b97f4a7c15ULL;
    return static_cast<double>((hash >> 37) & 0xffffu) / 65535.0;
}

size_t
mesh_physical_bucket(size_t mesh)
{
    const uint64_t hash = static_cast<uint64_t>(mesh) *
	11400714819323198485ULL + 0x9e3779b97f4a7c15ULL;
    return static_cast<size_t>((hash >> 17) % 10000);
}

double
mesh_physical_span(size_t mesh, double vehicleWidth)
{
    /* Deliberately decorrelate physical size from face count and database
     * order.  Real vehicle scenes contain a few hull-scale skins, hundreds of
     * appreciable structural pieces, ordinary components, and a long tail of
     * tiny fasteners.  A uniform 45..190-unit fixture inside a 10k-unit scene
     * hid priority inversions because every leaf occupied roughly the same
     * small number of pixels. */
    const size_t bucket = mesh_physical_bucket(mesh);
    const double unit = mesh_physical_unit(mesh);
    if (bucket < 10)
	return vehicleWidth * (0.35 + unit * 0.20); /* primary skins */
    if (bucket < 200)
	return vehicleWidth * (0.045 + unit * 0.14); /* structure */
    if (bucket < 1500)
	return 120.0 + unit * 680.0; /* subsystems */
    if (bucket < 6500)
	return 12.0 + unit * 190.0; /* ordinary components */
    return 0.8 + unit * 18.0; /* bolts, clips, and small fittings */
}

enum class MeshTopology {
    ClosedManifold,
    OpenSurface,
    NonManifold
};

MeshTopology
mesh_topology(size_t mesh)
{
    /* Vehicle-style workloads should overwhelmingly exercise the production
     * path: a closed, consistently oriented, shared-index surface.  Keep a
     * deterministic 5% adverse population so open and genuinely non-manifold
     * behavior does not regress without allowing it to dominate the visual or
     * performance result. */
    /* Permute each 40-leaf block independently.  This preserves the exact
     * ratio for full blocks without correlating topology with the repeating
     * face-count distribution. */
    const size_t block = mesh / 40;
    const size_t slot = ((mesh % 40) * 17 + block * 13) % 40;
    if (slot == 38)
	return MeshTopology::OpenSurface;
    if (slot == 39)
	return MeshTopology::NonManifold;
    return MeshTopology::ClosedManifold;
}

void
append_outward_face(std::vector<int> &faces,
    const std::vector<fastf_t> &vertices, int a, int b, int c)
{
    const fastf_t *av = &vertices[static_cast<size_t>(a) * 3];
    const fastf_t *bv = &vertices[static_cast<size_t>(b) * 3];
    const fastf_t *cv = &vertices[static_cast<size_t>(c) * 3];
    vect_t ab = {bv[X] - av[X], bv[Y] - av[Y], bv[Z] - av[Z]};
    vect_t ac = {cv[X] - av[X], cv[Y] - av[Y], cv[Z] - av[Z]};
    vect_t normal;
    VCROSS(normal, ab, ac);
    const point_t centroid = {
	(av[X] + bv[X] + cv[X]) / 3.0,
	(av[Y] + bv[Y] + cv[Y]) / 3.0,
	(av[Z] + bv[Z] + cv[Z]) / 3.0
    };
    if (VDOT(normal, centroid) < 0.0)
	std::swap(b, c);
    faces.push_back(a);
    faces.push_back(b);
    faces.push_back(c);
}

void
generate_closed_mesh(size_t mesh, size_t resolution,
    double spanX, double spanY, double spanZ,
    bool makeNonManifold, std::vector<fastf_t> &vertices,
    std::vector<int> &faces)
{
    /* A latitude/longitude ellipsoid gives us a closed shared-index surface
     * with approximately 2*resolution^2 faces: the same useful 500..100,000
     * distribution as the old heightfield fixture.  Deterministic ripples
     * make every asset's vertex data genuinely distinct. */
    const size_t longitudeSegments = resolution * 2;
    const size_t latitudeSegments = resolution / 2 + 1;
    const size_t interiorRings = latitudeSegments - 1;
    const size_t vertexCount =
	2 + interiorRings * longitudeSegments;
    vertices.assign(vertexCount * 3, 0.0);
    faces.clear();
    faces.reserve((2 * longitudeSegments * interiorRings +
	(makeNonManifold ? 1 : 0)) * 3);

    const double pi = std::acos(-1.0);
    const double phase = static_cast<double>(mesh) * 0.61803398875;
    const double radiusX = spanX * 0.5;
    const double radiusY = spanY * 0.5;
    const double radiusZ = spanZ * 0.5;
    vertices[X] = radiusX;
    const size_t southPole = vertexCount - 1;
    vertices[southPole * 3 + X] = -radiusX;

    for (size_t latitude = 1; latitude < latitudeSegments; ++latitude) {
	const double theta =
	    pi * static_cast<double>(latitude) /
	    static_cast<double>(latitudeSegments);
	const double sinTheta = std::sin(theta);
	const double cosTheta = std::cos(theta);
	for (size_t longitude = 0;
	     longitude < longitudeSegments; ++longitude) {
	    const double phi =
		2.0 * pi * static_cast<double>(longitude) /
		static_cast<double>(longitudeSegments);
	    const double ripple =
		1.0 +
		0.035 * std::sin(
		    phi * static_cast<double>(3 + mesh % 7) + phase) *
		    sinTheta * sinTheta +
		0.018 * std::cos(
		    theta * static_cast<double>(2 + mesh % 5) - phase);
	    const size_t vertex =
		1 + (latitude - 1) * longitudeSegments + longitude;
	    vertices[vertex * 3 + X] = radiusX * cosTheta * ripple;
	    vertices[vertex * 3 + Y] =
		radiusY * sinTheta * std::cos(phi) * ripple;
	    vertices[vertex * 3 + Z] =
		radiusZ * sinTheta * std::sin(phi) * ripple;
	}
    }

    for (size_t longitude = 0;
	 longitude < longitudeSegments; ++longitude) {
	const size_t next = (longitude + 1) % longitudeSegments;
	append_outward_face(faces, vertices, 0,
	    static_cast<int>(1 + longitude),
	    static_cast<int>(1 + next));
    }
    for (size_t ring = 0; ring + 1 < interiorRings; ++ring) {
	const size_t upper = 1 + ring * longitudeSegments;
	const size_t lower = upper + longitudeSegments;
	for (size_t longitude = 0;
	     longitude < longitudeSegments; ++longitude) {
	    const size_t next = (longitude + 1) % longitudeSegments;
	    append_outward_face(faces, vertices,
		static_cast<int>(upper + longitude),
		static_cast<int>(lower + longitude),
		static_cast<int>(lower + next));
	    append_outward_face(faces, vertices,
		static_cast<int>(upper + longitude),
		static_cast<int>(lower + next),
		static_cast<int>(upper + next));
	}
    }
    const size_t lastRing =
	1 + (interiorRings - 1) * longitudeSegments;
    for (size_t longitude = 0;
	 longitude < longitudeSegments; ++longitude) {
	const size_t next = (longitude + 1) % longitudeSegments;
	append_outward_face(faces, vertices,
	    static_cast<int>(lastRing + longitude),
	    static_cast<int>(southPole),
	    static_cast<int>(lastRing + next));
    }

    if (makeNonManifold && faces.size() >= 3) {
	/* A duplicate triangle gives each of its three edges three incident
	 * faces.  Positions and orientation remain plausible, but indexed
	 * topology analysis must reject this mesh as non-manifold. */
	faces.push_back(faces[0]);
	faces.push_back(faces[1]);
	faces.push_back(faces[2]);
    }
}

void
generate_open_mesh(size_t mesh, size_t gridCells,
    double spanX, double spanY, std::vector<fastf_t> &vertices,
    std::vector<int> &faces)
{
    const size_t side = gridCells + 1;
    const size_t vertexCount = side * side;
    const size_t faceCount = gridCells * gridCells * 2;
    vertices.resize(vertexCount * 3);
    faces.resize(faceCount * 3);
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
    const double amplitude = std::max(
	0.02, std::min(spanX, spanY) * 0.035);
    const size_t featureX = (mesh * 17 + 3) % side;
    const size_t featureY = (mesh * 29 + 5) % side;
    for (size_t y = 0; y < side; ++y) {
	for (size_t x = 0; x < side; ++x) {
	    const size_t vertex = y * side + x;
	    const double nx =
		static_cast<double>(x) / static_cast<double>(gridCells);
	    const double ny =
		static_cast<double>(y) / static_cast<double>(gridCells);
	    const double dx =
		static_cast<double>(x) - static_cast<double>(featureX);
	    const double dy =
		static_cast<double>(y) - static_cast<double>(featureY);
	    const double bump = std::exp(
		-(dx * dx + dy * dy) /
		std::max(4.0, static_cast<double>(gridCells) * 0.35));
	    vertices[vertex * 3 + X] = (nx - 0.5) * spanX;
	    vertices[vertex * 3 + Y] = (ny - 0.5) * spanY;
	    vertices[vertex * 3 + Z] =
		amplitude * std::sin(nx * 11.0 + phase) *
		    std::cos(ny * 9.0 - phase * 0.7) +
		bump * amplitude * 1.7;
	}
    }
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
    bu_setprogname(argv[0]);
    if (argc < 2 || argc > 4) {
	std::fprintf(stderr,
	    "Usage: %s output.g [mesh-count=5000] [max-grid-cells=224]\n",
	    argv[0]);
	return 2;
    }

    size_t meshCount = 5000;
    size_t maxGridCells = 224;
    if ((argc > 2 && !parse_positive(argv[2], 250000, meshCount)) ||
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
    const double vehicleWidth =
	static_cast<double>(std::max<size_t>(2, columns - 1)) * spacing;
    const double vehicleBeam = vehicleWidth * 0.32;
    const double vehicleHeight = vehicleWidth * 0.22;
    bool ok = true;
    size_t aggregateFaces = 0;
    size_t minimumFaces = std::numeric_limits<size_t>::max();
    size_t maximumFaces = 0;
    size_t closedManifoldCount = 0;
    size_t openSurfaceCount = 0;
    size_t nonManifoldCount = 0;

    for (size_t mesh = 0; mesh < meshCount && ok; ++mesh) {
	const size_t gridCells = mesh_grid_cells(mesh, maxGridCells);
	const double phase = static_cast<double>(mesh) * 0.61803398875;
	const double meshSpan = mesh_physical_span(mesh, vehicleWidth);
	const size_t physicalBucket = mesh_physical_bucket(mesh);
	const bool primarySkin = physicalBucket < 10;
	const bool structuralPanel =
	    physicalBucket >= 10 && physicalBucket < 200;
	const double profileUnit = mesh_physical_unit(mesh);
	const double meshSpanX = primarySkin ?
	    vehicleWidth * (0.22 + profileUnit * 0.14) : meshSpan;
	const double meshSpanY = primarySkin ?
	    vehicleBeam * (0.65 + profileUnit * 0.25) :
	    (structuralPanel ?
		meshSpan * (0.28 + profileUnit * 0.42) : meshSpan);
	const double meshSpanZ = primarySkin ?
	    vehicleHeight * (0.62 + profileUnit * 0.25) :
	    (structuralPanel ?
		meshSpan * (0.12 + profileUnit * 0.25) :
		meshSpan * (0.32 + profileUnit * 0.52));
	const MeshTopology topology = mesh_topology(mesh);
	std::vector<fastf_t> vertices;
	std::vector<int> faces;
	if (topology == MeshTopology::OpenSurface) {
	    generate_open_mesh(mesh, gridCells,
		meshSpanX, meshSpanY, vertices, faces);
	    ++openSurfaceCount;
	} else {
	    generate_closed_mesh(mesh, gridCells,
		meshSpanX, meshSpanY, meshSpanZ,
		topology == MeshTopology::NonManifold, vertices, faces);
	    if (topology == MeshTopology::ClosedManifold)
		++closedManifoldCount;
	    else
		++nonManifoldCount;
	}
	const size_t vertexCount = vertices.size() / 3;
	const size_t faceCount = faces.size() / 3;

	char leafName[64] = {0};
	const char *topologyName =
	    topology == MeshTopology::ClosedManifold ? "closed" :
	    (topology == MeshTopology::OpenSurface ? "open" :
		"nonmanifold");
	std::snprintf(leafName, sizeof(leafName), "unique_%s_%06zu.bot",
	    topologyName, mesh);
	const int botMode =
	    topology == MeshTopology::ClosedManifold ?
	    RT_BOT_SOLID : RT_BOT_SURFACE;
	if (mk_bot(wdbp, leafName, botMode, RT_BOT_CCW, 0,
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
	    rows > 1 ?
	    (static_cast<double>(row) -
		(static_cast<double>(rows) - 1.0) * 0.5) *
		(vehicleBeam / static_cast<double>(rows - 1)) : 0.0;
	const double tz =
	    layers > 1 ?
	    (static_cast<double>(layer) -
		(static_cast<double>(layers) - 1.0) * 0.5) *
		(vehicleHeight / static_cast<double>(layers - 1)) : 0.0;
	/* Hull-scale skins form a recognizable elongated envelope rather than
	 * inheriting the randomized panel placement.  This keeps visual failures
	 * diagnosable while the smaller population still stresses arbitrary
	 * placement, hierarchy, and visibility. */
	const double hullX = (profileUnit - 0.5) * vehicleWidth * 0.48;
	const double placedX = primarySkin ? hullX : tx;
	const double placedY = primarySkin ? 0.0 : ty;
	const double placedZ = primarySkin ? 0.0 :
	    (tz + 6.0 * std::sin(static_cast<double>(column) * 0.47) +
		4.0 * std::cos(static_cast<double>(row) * 0.53));
	/* Exercise projection, PCA, lighting, and culling with panels spanning
	 * all principal orientations plus a deterministic twist. */
	double twist = std::fmod(
	    static_cast<double>(mesh) * 137.50776405, 360.0);
	double rotateX = 7.0 * std::sin(phase);
	double rotateY = 7.0 * std::cos(phase * 0.7);
	if (primarySkin) {
	    twist = 0.0;
	    rotateX = 0.0;
	    rotateY = 0.0;
	} else if (structuralPanel) {
	    /* Structural members are mostly axis-aligned with a small authored
	     * cant, as opposed to the fully arbitrary smaller components. */
	    twist = static_cast<double>((mesh % 3) * 90);
	    rotateX = 3.0 * std::sin(phase);
	    rotateY = 3.0 * std::cos(phase * 0.7);
	} else {
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
	}
	bn_mat_angles(placement, rotateX, rotateY, twist);
	MAT_DELTAS(placement, placedX, placedY, placedZ);
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
    std::printf("generated %zu distinct meshes "
	"(%zu closed manifold, %zu open, %zu non-manifold), "
	"%zu..%zu faces each, %zu aggregate faces in %s\n",
	meshCount, closedManifoldCount, openSurfaceCount, nonManifoldCount,
	minimumFaces == std::numeric_limits<size_t>::max() ?
	    0 : minimumFaces,
	maximumFaces, aggregateFaces, argv[1]);
    return 0;
}
