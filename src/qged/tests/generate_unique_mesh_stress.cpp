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
#include <array>
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

enum class MeshShape : size_t {
    Body,
    Cylinder,
    Boom,
    Box,
    Panel,
    Dish,
    Irregular,
    Count
};

const char *
mesh_shape_name(MeshShape shape)
{
    switch (shape) {
	case MeshShape::Body: return "body";
	case MeshShape::Cylinder: return "cylinder";
	case MeshShape::Boom: return "boom";
	case MeshShape::Box: return "box";
	case MeshShape::Panel: return "panel";
	case MeshShape::Dish: return "dish";
	case MeshShape::Irregular: return "irregular";
	case MeshShape::Count: break;
    }
    return "unknown";
}

double
mesh_profile_unit(size_t mesh, uint64_t salt)
{
    uint64_t value = static_cast<uint64_t>(mesh) + salt;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    value ^= value >> 31;
    return static_cast<double>((value >> 11) & 0x1fffffffffffffULL) /
	static_cast<double>(0x1fffffffffffffULL);
}

bool
mesh_has_internal_orientation(size_t mesh)
{
    /* Combination transforms exercise instancing but leave each asset-local
     * AABB as tight as its OBB.  Internally orient a deterministic fifth of
     * the assets so scale fixtures also exercise cached PCA proxies found in
     * real databases whose authored local frames do not follow the part. */
    static const size_t populationPeriod = 5;
    static const size_t selectedResidue = 2;
    return mesh % populationPeriod == selectedResidue;
}

void
orient_mesh_vertices(size_t mesh, std::vector<fastf_t> &vertices)
{
    if (!mesh_has_internal_orientation(mesh))
	return;

    static const double minimumAngleDegrees = 17.0;
    static const double angleRangeDegrees = 31.0;
    mat_t rotation;
    bn_mat_angles(rotation,
	minimumAngleDegrees + angleRangeDegrees *
	    mesh_profile_unit(mesh, 0xa4093822299f31d0ULL),
	minimumAngleDegrees + angleRangeDegrees *
	    mesh_profile_unit(mesh, 0x082efa98ec4e6c89ULL),
	minimumAngleDegrees + angleRangeDegrees *
	    mesh_profile_unit(mesh, 0x452821e638d01377ULL));
    for (size_t offset = 0; offset + 2 < vertices.size(); offset += 3) {
	point_t source;
	point_t oriented;
	VSET(source, vertices[offset], vertices[offset + 1],
	    vertices[offset + 2]);
	MAT4X3PNT(oriented, rotation, source);
	vertices[offset] = oriented[X];
	vertices[offset + 1] = oriented[Y];
	vertices[offset + 2] = oriented[Z];
    }
}

MeshShape
mesh_shape(size_t mesh, MeshTopology topology, size_t physicalBucket)
{
    /* Keep related leaves related.  Each 256-leaf run occupies one of twelve
     * recognizable spacecraft/vehicle assemblies, so the balanced database
     * hierarchy contains branches with very different visual significance.
     * A uniformly shuffled fixture can make a prefix-biased scheduler look
     * representative by accident. */
    const size_t assembly = (mesh / 256) % 12;
    const size_t selector = static_cast<size_t>(
	mesh_profile_unit(mesh, 0x243f6a8885a308d3ULL) * 100.0);

    if (topology == MeshTopology::OpenSurface)
	return (assembly == 6 || assembly == 7 || selector < 28) ?
	    MeshShape::Dish : MeshShape::Panel;

    /* Hull-scale parts define the overall silhouette.  Their family follows
     * their assembly instead of degenerating into a collection of uniformly
     * scaled ellipsoids. */
    if (physicalBucket < 10) {
	if (assembly == 2 || assembly == 3)
	    return MeshShape::Panel;
	if (assembly == 4 || assembly == 5)
	    return MeshShape::Boom;
	return assembly == 6 ? MeshShape::Cylinder : MeshShape::Body;
    }

    if (assembly == 2 || assembly == 3) {
	if (selector < 62) return MeshShape::Panel;
	if (selector < 82) return MeshShape::Boom;
	if (selector < 94) return MeshShape::Box;
	return MeshShape::Cylinder;
    }
    if (assembly == 4 || assembly == 5) {
	if (selector < 58) return MeshShape::Boom;
	if (selector < 80) return MeshShape::Cylinder;
	if (selector < 93) return MeshShape::Box;
	return MeshShape::Irregular;
    }
    if (assembly == 6 || assembly == 7) {
	if (selector < 30) return MeshShape::Cylinder;
	if (selector < 52) return MeshShape::Body;
	if (selector < 75) return MeshShape::Box;
	if (selector < 88) return MeshShape::Boom;
	return MeshShape::Irregular;
    }
    if (selector < 24) return MeshShape::Box;
    if (selector < 43) return MeshShape::Cylinder;
    if (selector < 60) return MeshShape::Boom;
    if (selector < 74) return MeshShape::Panel;
    if (selector < 92) return MeshShape::Body;
    return MeshShape::Irregular;
}

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
    MeshShape shape, bool makeNonManifold, std::vector<fastf_t> &vertices,
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
	    const double sphereX = cosTheta;
	    const double sphereY = sinTheta * std::cos(phi);
	    const double sphereZ = sinTheta * std::sin(phi);
	    const auto signedPower = [](double value, double exponent) {
		return std::copysign(std::pow(std::fabs(value), exponent), value);
	    };
	    double shapeX = sphereX;
	    double shapeY = sphereY;
	    double shapeZ = sphereZ;
	    switch (shape) {
		case MeshShape::Cylinder:
		case MeshShape::Boom:
		    /* A low axial exponent produces a long cylindrical middle
		     * and rounded caps without changing the closed sphere
		     * topology or introducing coincident seam vertices. */
		    shapeX = signedPower(sphereX,
			shape == MeshShape::Boom ? 0.18 : 0.30);
		    break;
		case MeshShape::Box:
		case MeshShape::Panel:
		    shapeX = signedPower(sphereX, 0.24);
		    shapeY = signedPower(sphereY, 0.24);
		    shapeZ = signedPower(sphereZ, 0.24);
		    break;
		case MeshShape::Body:
		    shapeX = signedPower(sphereX, 0.58);
		    break;
		case MeshShape::Irregular: {
		    const double skew = 0.10 * sphereX * sphereY;
		    shapeY += skew + 0.045 * std::sin(3.0 * phi + phase) *
			sinTheta * sinTheta;
		    shapeZ += 0.07 * sphereX * sphereY;
		    break;
		}
		case MeshShape::Dish:
		case MeshShape::Count:
		    break;
	    }
	    vertices[vertex * 3 + X] = radiusX * shapeX * ripple;
	    vertices[vertex * 3 + Y] = radiusY * shapeY * ripple;
	    vertices[vertex * 3 + Z] = radiusZ * shapeZ * ripple;
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
generate_open_dish(size_t mesh, size_t resolution,
    double spanX, double spanY, std::vector<fastf_t> &vertices,
    std::vector<int> &faces)
{
    const size_t radialSegments = resolution / 2 + 1;
    const size_t angularSegments = resolution * 2;
    const size_t vertexCount = 1 + radialSegments * angularSegments;
    vertices.assign(vertexCount * 3, 0.0);
    faces.clear();
    faces.reserve((angularSegments +
	(radialSegments - 1) * angularSegments * 2) * 3);

    const double pi = std::acos(-1.0);
    const double phase = static_cast<double>(mesh) * 0.61803398875;
    const double depth = std::min(spanX, spanY) *
	(0.10 + 0.08 * mesh_profile_unit(mesh, 0xa4093822299f31d0ULL));
    for (size_t radial = 1; radial <= radialSegments; ++radial) {
	const double radius = static_cast<double>(radial) /
	    static_cast<double>(radialSegments);
	for (size_t angular = 0; angular < angularSegments; ++angular) {
	    const double phi = 2.0 * pi * static_cast<double>(angular) /
		static_cast<double>(angularSegments);
	    const size_t vertex =
		1 + (radial - 1) * angularSegments + angular;
	    vertices[vertex * 3 + X] = spanX * 0.5 * radius * std::cos(phi);
	    vertices[vertex * 3 + Y] = spanY * 0.5 * radius * std::sin(phi);
	    vertices[vertex * 3 + Z] = depth * radius * radius +
		depth * 0.025 * std::sin(5.0 * phi + phase) * radius;
	}
    }

    for (size_t angular = 0; angular < angularSegments; ++angular) {
	const size_t next = (angular + 1) % angularSegments;
	faces.push_back(0);
	faces.push_back(static_cast<int>(1 + angular));
	faces.push_back(static_cast<int>(1 + next));
    }
    for (size_t radial = 1; radial < radialSegments; ++radial) {
	const size_t inner = 1 + (radial - 1) * angularSegments;
	const size_t outer = inner + angularSegments;
	for (size_t angular = 0; angular < angularSegments; ++angular) {
	    const size_t next = (angular + 1) % angularSegments;
	    faces.push_back(static_cast<int>(inner + angular));
	    faces.push_back(static_cast<int>(outer + angular));
	    faces.push_back(static_cast<int>(outer + next));
	    faces.push_back(static_cast<int>(inner + angular));
	    faces.push_back(static_cast<int>(outer + next));
	    faces.push_back(static_cast<int>(inner + next));
	}
    }
}

void
generate_open_mesh(size_t mesh, size_t gridCells,
    double spanX, double spanY, MeshShape shape,
    std::vector<fastf_t> &vertices,
    std::vector<int> &faces)
{
    if (shape == MeshShape::Dish) {
	generate_open_dish(mesh, gridCells, spanX, spanY, vertices, faces);
	return;
    }
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
    size_t internallyOrientedCount = 0;
    std::array<size_t, static_cast<size_t>(MeshShape::Count)> shapeCounts = {};
    double minimumAspectRatio = 1.0;
    double minimumPhysicalSpan = std::numeric_limits<double>::max();
    double maximumPhysicalSpan = 0.0;

    for (size_t mesh = 0; mesh < meshCount && ok; ++mesh) {
	const size_t gridCells = mesh_grid_cells(mesh, maxGridCells);
	const double phase = static_cast<double>(mesh) * 0.61803398875;
	const double meshSpan = mesh_physical_span(mesh, vehicleWidth);
	const size_t physicalBucket = mesh_physical_bucket(mesh);
	const bool primarySkin = physicalBucket < 10;
	const bool structuralPanel =
	    physicalBucket >= 10 && physicalBucket < 200;
	const double profileUnit = mesh_physical_unit(mesh);
	const MeshTopology topology = mesh_topology(mesh);
	MeshShape shape = mesh_shape(mesh, topology, physicalBucket);
	/* Preserve the established 5k smooth-zoom target as a detailed body;
	 * its name, source density, placement, and useful center are harness
	 * invariants rather than part of the statistical profile. */
	if (mesh == 199)
	    shape = MeshShape::Body;
	const double aspectUnit =
	    mesh_profile_unit(mesh, 0x13198a2e03707344ULL);
	double meshSpanX = meshSpan;
	double meshSpanY = meshSpan * (0.48 + profileUnit * 0.42);
	double meshSpanZ = meshSpan * (0.36 + aspectUnit * 0.45);
	switch (shape) {
	    case MeshShape::Body:
		if (primarySkin) {
		    meshSpanX = vehicleWidth * (0.22 + profileUnit * 0.14);
		    meshSpanY = vehicleBeam * (0.65 + aspectUnit * 0.25);
		    meshSpanZ = vehicleHeight * (0.62 + profileUnit * 0.25);
		}
		break;
	    case MeshShape::Cylinder:
		meshSpanY = meshSpan * (0.18 + profileUnit * 0.30);
		meshSpanZ = meshSpanY * (0.82 + aspectUnit * 0.30);
		break;
	    case MeshShape::Boom:
		meshSpanY = meshSpan * (0.018 + profileUnit * 0.065);
		meshSpanZ = meshSpan * (0.018 + aspectUnit * 0.065);
		break;
	    case MeshShape::Box:
		meshSpanY = meshSpan * (0.30 + profileUnit * 0.55);
		meshSpanZ = meshSpan * (0.14 + aspectUnit * 0.46);
		break;
	    case MeshShape::Panel:
		meshSpanY = meshSpan * (0.42 + profileUnit * 0.70);
		meshSpanZ = meshSpan * (0.006 + aspectUnit * 0.026);
		break;
	    case MeshShape::Dish:
		meshSpanY = meshSpan * (0.68 + profileUnit * 0.34);
		meshSpanZ = meshSpan * (0.10 + aspectUnit * 0.08);
		break;
	    case MeshShape::Irregular:
		meshSpanY = meshSpan * (0.38 + profileUnit * 0.70);
		meshSpanZ = meshSpan * (0.28 + aspectUnit * 0.62);
		break;
	    case MeshShape::Count:
		break;
	}
	if (structuralPanel && shape != MeshShape::Panel &&
		shape != MeshShape::Boom) {
	    meshSpanY *= 0.72;
	    meshSpanZ *= 0.62;
	}
	if (mesh == 199) {
	    meshSpanX = meshSpan;
	    meshSpanY = meshSpan;
	    meshSpanZ = meshSpan * (0.32 + profileUnit * 0.52);
	}
	++shapeCounts[static_cast<size_t>(shape)];
	const double shortestSpan = std::min(meshSpanX,
	    std::min(meshSpanY, meshSpanZ));
	const double longestSpan = std::max(meshSpanX,
	    std::max(meshSpanY, meshSpanZ));
	minimumAspectRatio = std::min(minimumAspectRatio,
	    shortestSpan / std::max(0.001, longestSpan));
	minimumPhysicalSpan = std::min(minimumPhysicalSpan, longestSpan);
	maximumPhysicalSpan = std::max(maximumPhysicalSpan, longestSpan);
	std::vector<fastf_t> vertices;
	std::vector<int> faces;
	if (topology == MeshTopology::OpenSurface) {
	    generate_open_mesh(mesh, gridCells,
		meshSpanX, meshSpanY, shape, vertices, faces);
	    ++openSurfaceCount;
	} else {
	    generate_closed_mesh(mesh, gridCells,
		meshSpanX, meshSpanY, meshSpanZ, shape,
		topology == MeshTopology::NonManifold, vertices, faces);
	    if (topology == MeshTopology::ClosedManifold)
		++closedManifoldCount;
	    else
		++nonManifoldCount;
	}
	orient_mesh_vertices(mesh, vertices);
	if (mesh_has_internal_orientation(mesh))
	    ++internallyOrientedCount;
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
	double placedX = primarySkin ? hullX : tx;
	double placedY = primarySkin ? 0.0 : ty;
	double placedZ = primarySkin ? 0.0 :
	    (tz + 6.0 * std::sin(static_cast<double>(column) * 0.47) +
		4.0 * std::cos(static_cast<double>(row) * 0.53));
	if (mesh != 199 && !primarySkin) {
	    /* A dozen authored macro-assemblies replace the old perturbed lattice.
	     * Consecutive hierarchy branches are spatially and visually coherent,
	     * while deterministic local coordinates keep all assets distinct. */
	    const size_t assembly = (mesh / 256) % 12;
	    const double localX =
		mesh_profile_unit(mesh, 0x082efa98ec4e6c89ULL) * 2.0 - 1.0;
	    const double localY =
		mesh_profile_unit(mesh, 0x452821e638d01377ULL) * 2.0 - 1.0;
	    const double localZ =
		mesh_profile_unit(mesh, 0xbe5466cf34e90c6cULL) * 2.0 - 1.0;
	    switch (assembly) {
		case 0:
		case 1: {
		    const double shellAngle = localY * std::acos(-1.0);
		    placedX = localX * vehicleWidth * 0.34;
		    placedY = std::cos(shellAngle) * vehicleBeam * 0.26;
		    placedZ = std::sin(shellAngle) * vehicleHeight * 0.30;
		    break;
		}
		case 2:
		case 3: {
		    const double side = assembly == 2 ? -1.0 : 1.0;
		    placedX = localX * vehicleWidth * 0.15;
		    placedY = side * vehicleWidth * (0.22 + 0.28 *
			(0.5 + 0.5 * localY));
		    placedZ = localZ * vehicleHeight * 0.10;
		    break;
		}
		case 4:
		case 5: {
		    const double side = assembly == 4 ? -1.0 : 1.0;
		    const double along = 0.5 + 0.5 * localX;
		    placedX = side * along * vehicleWidth * 0.48;
		    placedY = localY * vehicleBeam * (0.08 + 0.20 * along);
		    placedZ = localZ * vehicleHeight * (0.08 + 0.18 * along);
		    break;
		}
		case 6:
		case 7:
		    placedX = (assembly == 6 ? -0.18 : 0.18) * vehicleWidth +
			localX * vehicleWidth * 0.09;
		    placedY = localY * vehicleBeam * 0.30;
		    placedZ = vehicleHeight * 0.28 +
			localZ * vehicleHeight * 0.18;
		    break;
		case 8:
		case 9:
		    placedX = (assembly == 8 ? -0.31 : 0.31) * vehicleWidth +
			localX * vehicleWidth * 0.11;
		    placedY = localY * vehicleBeam * 0.34;
		    placedZ = -vehicleHeight * 0.20 +
			localZ * vehicleHeight * 0.20;
		    break;
		default:
		    placedX = localX * vehicleWidth * 0.38;
		    placedY = localY * vehicleBeam * 0.42;
		    placedZ = localZ * vehicleHeight * 0.45;
		    break;
	    }
	}
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
	} else if (structuralPanel || shape == MeshShape::Panel ||
		shape == MeshShape::Boom) {
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
    if (meshCount >= 1000) {
	bool profileOk = minimumAspectRatio < 0.05 &&
	    maximumPhysicalSpan > minimumPhysicalSpan * 100.0;
	profileOk = profileOk && internallyOrientedCount >= meshCount / 6;
	const size_t minimumShapePopulation = meshCount / 200;
	for (size_t shape = 0;
	     shape < static_cast<size_t>(MeshShape::Count); ++shape) {
	    profileOk = profileOk &&
		shapeCounts[shape] >= minimumShapePopulation &&
		shapeCounts[shape] < meshCount / 2;
	}
	if (!profileOk) {
	    std::fprintf(stderr,
		"generated profile lacks required size/aspect/shape/orientation diversity\n");
	    return 1;
	}
    }
    std::printf("generated %zu distinct meshes "
	"(%zu closed manifold, %zu open, %zu non-manifold), "
	"%zu..%zu faces each, %zu aggregate faces in %s\n",
	meshCount, closedManifoldCount, openSurfaceCount, nonManifoldCount,
	minimumFaces == std::numeric_limits<size_t>::max() ?
	    0 : minimumFaces,
	maximumFaces, aggregateFaces, argv[1]);
    std::printf("shape profile:");
    for (size_t shape = 0;
	 shape < static_cast<size_t>(MeshShape::Count); ++shape)
	std::printf(" %s=%zu", mesh_shape_name(static_cast<MeshShape>(shape)),
	    shapeCounts[shape]);
    std::printf(", internally-oriented=%zu, aspect-min=%.4f, "
	"physical-span-ratio=%.1f\n", internallyOrientedCount,
	minimumAspectRatio,
	maximumPhysicalSpan / std::max(0.001, minimumPhysicalSpan));
    return 0;
}
