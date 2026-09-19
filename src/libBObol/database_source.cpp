/*                D A T A B A S E _ S O U R C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"
#include "scalar_publication_private.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BDrawCache.h"
#include "BObol/BEvaluatedPoints.h"
#include "BObol/BExportAction.h"
#include "BObol/BLodMeshShape.h"
#include "BObol/BLodRealization.h"
#include "BObol/BLodService.h"
#include "BObol/BMaterialObject.h"
#include "BObol/BMeasureAction.h"
#include "BObol/BMeshLodCache.h"
#include "BObol/BMeshShape.h"
#include "BObol/BPickDetail.h"
#include "BObol/BSnapAction.h"
#include "BObol/BViewLod.h"
#include "BObol/BViewQuery.h"
#include "BObol/BVListShape.h"
#include "cad_assembly_private.h"
#include "cad_normals_private.h"
#include "cad_publication_private.h"
#include "compact_occurrence_registry_private.h"

#include "database_source_private.h"
#include "database_source_mesh_geometry_private.h"
#include "database_source_presentation_private.h"
#include "database_source_realization.h"
#include "identity_counter_private.h"
#include "performance_private.h"
#include "serialized_bot_source_private.h"
#include "transaction_fault_private.h"

#include "bg/line_layer.h"
#include "bg/pca.h"
#include "bg/trimesh.h"
#include "bg/vlist.h"
#include "bu/app.h"
#include "bu/color.h"
#include "bu/cv.h"
#include "bu/file.h"
#include "bu/hash.h"
#include "bu/list.h"
#include "bu/mapped_file.h"
#include "bu/parallel.h"
#include "bu/str.h"
#include "bu/datetime.h"
#include "bu/vls.h"
#include "nmg.h"
#include "raytrace.h"
#include "rt/func.h"
#include "rt/global.h"
#include "rt/db4.h"
#include "rt/nongeom.h"
#include "rt/db_fullpath.h"
#include "rt/display_bounds.h"
#include "rt/eval_wireframe.h"
#include "rt/primitives/annot.h"
#include "rt/tree.h"
#include "rt/vlist.h"
#include "rt/view.h"
#include "wdb.h"

#include <Inventor/SbName.h>
#include <Inventor/fields/SoFieldData.h>
#include <Inventor/fields/SoMField.h>
#include <Inventor/fields/SoSField.h>
#include <Inventor/tools/SbModernUtils.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/nodes/SoGroup.h>

#include <Inventor/nodes/SoMatrixTransform.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/misc/SoNotification.h>
#include <Obol/cad/CadProjectedProxy.h>
#include <Inventor/sensors/SoFieldSensor.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <exception>
#include <initializer_list>
#include <Inventor/misc/SoChildList.h>
#include <inttypes.h>
#include <limits.h>
#include <limits>
#include <map>
#include <math.h>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <set>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

void realized_vlist_shape_summary(const SoBRLVListShape *shape,
	BObolRealizedShapeSummary &summary);
void realized_mesh_shape_summary(const SoBRLMeshShape *shape,
	BObolRealizedShapeSummary &summary);
BObolDatabaseSourcePublishState::BObolDatabaseSourcePublishState(void) :
    sourceInstanceKey(NULL),
    sourcePath(NULL),
    sourceRepresentationKey(NULL),
    targetGroupPath(NULL),
    database(NULL),
    drawMode(SoBRLDatabaseSource::WIREFRAME),
    representationMode(SoBRLDatabaseSource::REPRESENTATION_DEFAULT),
    sourceRevisionValid(FALSE),
    sourceRevision(0),
    inputsRevision(0),
    visible(TRUE),
    selected(FALSE),
    highlighted(FALSE),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0f),
    colorOverride(FALSE),
    color(1.0f, 1.0f, 1.0f),
    materialColorValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialRevision(0),
    materialPolicyValid(FALSE),
    materialPolicy(SoBRLDatabaseSource::MATERIAL_INHERIT),
    roleFlagsValid(FALSE),
    roleFlags(SoBRLDatabaseSource::REALIZATION_ROLE_NONE),
    realizationStateValid(FALSE),
    realizationStatus(SoBRLDatabaseSource::UNREALIZED),
    realizedSourceRevision(0),
    realizedInputsRevision(0),
    staleReason(SoBRLDatabaseSource::STALE_SOURCE),
    realizationDiagnostic(NULL),
    viewPolicyValid(FALSE),
    viewDependent(FALSE),
    csgLodEnabled(FALSE),
    meshLodEnabled(FALSE),
    viewScale(0.0f),
    lodScale(1.0f),
    viewWidth(0),
    viewHeight(0),
    botThreshold(0),
    curveScale(0.0f),
    pointScale(0.0f),
    placementValid(FALSE),
    drawMatrixValid(FALSE),
    drawMatrix(SbMatrix::identity()),
    drawCenterValid(FALSE),
    drawCenter(0.0f, 0.0f, 0.0f),
    drawSizeValid(FALSE),
    drawSize(0.0f)
{
}


static int cad_source_mesh_request_from_bot(
	BObolSourceMeshRequest &request, const struct rt_bot_internal *bot);
static uint64_t
database_source_encoded_source_bytes(const struct db_i *dbip,
	const struct directory *dp)
{
    if (!dp || dp->d_len <= 0)
	return 0;
    uint64_t encodedBytes = static_cast<uint64_t>(dp->d_len);
    if (dbip && db_version(dbip) < 5) {
	const uint64_t recordBytes = sizeof(union record);
	if (encodedBytes > UINT64_MAX / recordBytes)
	    return UINT64_MAX;
	encodedBytes *= recordBytes;
    }
    return encodedBytes;
}

const SbString &
compact_instance_identity(const BObolCompactInstanceEntry &entry)
{
    return entry.instanceKey;
}

SbBox3f
compact_part_geometry_bounds(
    const std::shared_ptr<const Obol::PartGeometry> &geometry)
{
    SbBox3f bounds;
    bounds.makeEmpty();
    if (!geometry)
	return bounds;
    if (geometry->conservativeBounds &&
	!geometry->conservativeBounds->isEmpty())
	bounds.extendBy(*geometry->conservativeBounds);
    if (geometry->points)
	bounds.extendBy(geometry->points->bounds);
    if (geometry->wire)
	bounds.extendBy(geometry->wire->bounds);
    if (geometry->shaded)
	bounds.extendBy(geometry->shaded->bounds);
    return bounds;
}

SO_NODE_SOURCE(SoBRLDatabaseSource);

static std::atomic<uint64_t> database_source_next_handle_id(1);

uint64_t
database_source_handle_id(void)
{
    return bobol_atomic_nonzero_identity_take(
	database_source_next_handle_id, std::memory_order_seq_cst,
	std::memory_order_seq_cst);
}

int
database_source_float_different(float a, float b)
{
    return fabsf(a - b) > 1.0e-6f;
}

int
database_source_color_equal(const SbColor &a, const SbColor &b)
{
    return !database_source_float_different(a[0], b[0]) &&
	   !database_source_float_different(a[1], b[1]) &&
	   !database_source_float_different(a[2], b[2]);
}

static int
database_source_vec3f_equal(const SbVec3f &a, const SbVec3f &b)
{
    return !database_source_float_different(a[0], b[0]) &&
	   !database_source_float_different(a[1], b[1]) &&
	   !database_source_float_different(a[2], b[2]);
}

int
database_source_string_equal(const SbString &a, const char *b)
{
    return bu_strcmp(a.getString(), b ? b : "") == 0;
}

template <typename FieldT>
static void
database_source_assign_string(FieldT &field, const char *value)
{
    const char *nextValue = value ? value : "";
    if (database_source_string_equal(field.getValue(), nextValue))
	return;

    const std::string stableValue(nextValue);
    field = stableValue.c_str();
}

template <typename FieldT>
static void
database_source_assign_string(FieldT &field, const SbString &value)
{
    database_source_assign_string(field, value.getString());
}

static SbBox3f
database_source_box_from_minmax(const SbVec3f &bmin, const SbVec3f &bmax)
{
    SbBox3f bounds;
    bounds.makeEmpty();
    bounds.extendBy(bmin);
    bounds.extendBy(bmax);
    return bounds;
}

SbBox3f
database_source_transform_bounds(const SbBox3f &bounds,
				 const SbMatrix &matrix)
{
    SbBox3f transformed;
    transformed.makeEmpty();
    if (bounds.isEmpty())
	return transformed;

    const SbVec3f bmin = bounds.getMin();
    const SbVec3f bmax = bounds.getMax();
    for (int xi = 0; xi < 2; xi++) {
	for (int yi = 0; yi < 2; yi++) {
	    for (int zi = 0; zi < 2; zi++) {
		const SbVec3f corner(
		    xi ? bmax[0] : bmin[0],
		    yi ? bmax[1] : bmin[1],
		    zi ? bmax[2] : bmin[2]);
		SbVec3f transformedCorner;
		matrix.multVecMatrix(corner, transformedCorner);
		transformed.extendBy(transformedCorner);
	    }
	}
    }

    return transformed;
}

static const char *
database_source_material_skip_leading_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

std::string
database_source_db_path_without_instance_suffixes(const char *path)
{
    std::string lookupPath;
    if (!path)
	return lookupPath;

    for (const char *cp = database_source_material_skip_leading_slash(path);
	 *cp;
	 cp++) {
	if (*cp == '@' && cp[1] >= '0' && cp[1] <= '9') {
	    while (cp[1] && cp[1] != '/')
		cp++;
	    continue;
	}
	lookupPath.push_back(*cp);
    }
    return lookupPath;
}

SbBool
bobol_database_source_path_material_color(
    struct db_i *dbip,
    const char *path,
    SbColor &color)
{
    if (!dbip || !path || !path[0])
	return FALSE;

    const std::string lookupPath =
	database_source_db_path_without_instance_suffixes(path);
    if (lookupPath.empty())
	return FALSE;

    struct db_full_path fullpath;
    db_full_path_init(&fullpath);
    if (db_string_to_path(&fullpath, dbip, lookupPath.c_str()) != 0) {
	db_free_full_path(&fullpath);
	return FALSE;
    }

    const SbBool matched =
	bobol_database_source_fullpath_material_color(dbip, &fullpath,
						       color);
    db_free_full_path(&fullpath);
    return matched;
}

namespace {

struct BObolMaterialObjectState {
    BObolMaterialObjectState(void) :
	isCombination(false),
	isAnnotation(false),
	rgbValid(false),
	inherit(DB_INH_LOWER),
	isRegion(false),
	regionId(-1),
	airCode(0),
	materialId(0),
	los(0)
    {
	rgb[0] = rgb[1] = rgb[2] = 0;
    }

    bool isCombination;
    bool isAnnotation;
    bool rgbValid;
    unsigned char rgb[3];
    int inherit;
    bool isRegion;
    int regionId;
    int airCode;
    int materialId;
    int los;
    std::string shader;
};

struct BObolMaterialPathState {
    BObolMaterialPathState(void) :
	explicitColorValid(false),
	isAnnotation(false),
	colorInherit(DB_INH_LOWER),
	inRegion(false),
	regionId(-1),
	airCode(0),
	materialId(0),
	los(0),
	color(1.0f, 0.0f, 0.0f)
    {
    }

    bool explicitColorValid;
    bool isAnnotation;
    int colorInherit;
    bool inRegion;
    int regionId;
    int airCode;
    int materialId;
    int los;
    std::string shader;
    SbColor color;
};

class BObolMaterialColorSweep {
public:
    explicit BObolMaterialColorSweep(struct db_i *database) : dbip(database)
    {
    }

    bool resolve(const char *sourcePath, BObolMaterialPathState &resolved)
    {
	if (!this->dbip || !sourcePath || !sourcePath[0])
	    return false;

	const std::string path =
	    database_source_db_path_without_instance_suffixes(sourcePath);
	const char *cp = path.c_str();
	while (*cp == '/')
	    cp++;
	if (!cp[0])
	    return false;

	BObolMaterialPathState state;
	std::string prefix;
	while (*cp) {
	    const char *slash = strchr(cp, '/');
	    const size_t length = slash ? static_cast<size_t>(slash - cp) :
		strlen(cp);
	    if (!length) {
		cp = slash ? slash + 1 : cp + length;
		continue;
	    }

	    const std::string component(cp, length);
	    if (!prefix.empty())
		prefix.push_back('/');
	    prefix.append(component);

	    std::unordered_map<std::string,
		BObolMaterialPathState>::const_iterator cached =
		this->pathStates.find(prefix);
	    if (cached != this->pathStates.end()) {
		state = cached->second;
	    } else {
		struct directory *dp = db_lookup(this->dbip,
			component.c_str(), LOOKUP_QUIET);
		if (!dp)
		    return false;
		this->applyObject(state, dp);
		this->pathStates.emplace(prefix, state);
	    }

	    if (!slash)
		break;
	    cp = slash + 1;
	}

	if (!state.explicitColorValid && !state.isAnnotation && state.regionId >= 0)
	    state.color = this->regionColor(state.regionId);
	resolved = state;
	return true;
    }

    bool resolve(const struct db_full_path *sourcePath,
	BObolMaterialPathState &resolved)
    {
	if (!this->dbip || !sourcePath || sourcePath->fp_len == 0 ||
	    !sourcePath->fp_names)
	    return false;

	BObolMaterialPathNode *node = &this->pathRoot;
	for (size_t i = 0; i < sourcePath->fp_len; ++i) {
	    struct directory *dp = sourcePath->fp_names[i];
	    if (!dp)
		return false;
	    std::unique_ptr<BObolMaterialPathNode> &child =
		node->children[dp];
	    if (!child) {
		child.reset(new BObolMaterialPathNode);
		child->state = node->state;
		this->applyObject(child->state, dp);
	    }
	    node = child.get();
	}

	resolved = node->state;
	if (!resolved.explicitColorValid && !resolved.isAnnotation && resolved.regionId >= 0)
	    resolved.color = this->regionColor(resolved.regionId);
	return true;
    }

private:
    struct BObolMaterialPathNode {
	BObolMaterialPathState state;
	std::unordered_map<struct directory *,
	    std::unique_ptr<BObolMaterialPathNode>> children;
    };

    const BObolMaterialObjectState &objectState(struct directory *dp)
    {
	std::unordered_map<struct directory *,
	    BObolMaterialObjectState>::iterator found =
	    this->objectStates.find(dp);
	if (found != this->objectStates.end())
	    return found->second;

	BObolMaterialObjectState state;
	state.isCombination = (dp->d_flags & RT_DIR_COMB) != 0;
	if (state.isCombination) {
	    struct rt_db_internal intern;
	    RT_DB_INTERNAL_INIT(&intern);
	    if (rt_db_get_internal(&intern, dp, this->dbip, NULL) >= 0) {
		const std::unique_ptr<rt_db_internal, decltype(&rt_db_free_internal)> imported(&intern, rt_db_free_internal);
		if (intern.idb_type == ID_COMBINATION && intern.idb_ptr) {
		    const struct rt_comb_internal *comb =
			static_cast<const struct rt_comb_internal *>(intern.idb_ptr);
		    RT_CK_COMB(comb);
		    state.rgbValid = comb->rgb_valid == 1;
		    if (state.rgbValid) {
			state.rgb[0] = comb->rgb[0];
			state.rgb[1] = comb->rgb[1];
			state.rgb[2] = comb->rgb[2];
		    }
		    state.inherit = comb->inherit;
		    state.isRegion = comb->region_flag != 0;
		    state.regionId = comb->region_id;
		    state.airCode = comb->aircode;
		    state.materialId = comb->GIFTmater;
		    state.los = comb->los;
		    state.shader = bu_vls_cstr(&comb->shader);
		}
	    }
	} else {
	    state.isAnnotation = dp->d_major_type == DB5_MAJORTYPE_BRLCAD &&
		dp->d_minor_type == DB5_MINORTYPE_BRLCAD_ANNOT;
	    struct bu_attribute_value_set attributes = BU_AVS_INIT_ZERO;
	    if (db5_get_attributes(this->dbip, &attributes, dp) == 0) {
		const char *attribute = bu_avs_get(&attributes,
		    db5_standard_attribute(ATTR_COLOR));
		int components[3];
		/* Match primitive wireframe and shaded color parsing, including
		 * clamping; malformed attributes retain the parent color. */
		if (attribute && sscanf(attribute, "%3i%*c%3i%*c%3i",
		    &components[0], &components[1], &components[2]) == 3 &&
		    components[0] >= 0 && components[1] >= 0 && components[2] >= 0) {
		    const int maxColorComponent = 255;
		    for (int i = 0; i < 3; ++i)
			state.rgb[i] = static_cast<unsigned char>(
			    std::min(components[i], maxColorComponent));
		    state.rgbValid = true;
		}
	    }
	    bu_avs_free(&attributes);
	}
	return this->objectStates.emplace(dp, state).first->second;
    }

    void applyObject(BObolMaterialPathState &pathState,
	struct directory *dp)
    {
	const BObolMaterialObjectState &comb = this->objectState(dp);
	pathState.isAnnotation = comb.isAnnotation;
	if (comb.isAnnotation && !pathState.explicitColorValid)
	    pathState.color = SbColor(1.0f, 1.0f, 1.0f);

	/* Primitive attributes override inherited combination colors in every
	 * drawing mode.  Combination inheritance and region boundaries still
	 * control colors specified by other combinations. */
	if (comb.rgbValid && (!comb.isCombination ||
	    (!pathState.inRegion && pathState.colorInherit == DB_INH_LOWER))) {
	    pathState.color = SbColor(
		static_cast<float>(comb.rgb[0]) / 255.0f,
		static_cast<float>(comb.rgb[1]) / 255.0f,
		static_cast<float>(comb.rgb[2]) / 255.0f);
	    pathState.explicitColorValid = true;
	    pathState.colorInherit = comb.inherit;
	}
	if (!pathState.inRegion && comb.isRegion) {
	    pathState.regionId = comb.regionId;
	    pathState.airCode = comb.airCode;
	    pathState.materialId = comb.materialId;
	    pathState.los = comb.los;
	    pathState.shader = comb.shader;
	    pathState.inRegion = true;
	}
    }

    SbColor regionColor(int regionId)
    {
	std::unordered_map<int, SbColor>::const_iterator found =
	    this->regionColors.find(regionId);
	if (found != this->regionColors.end())
	    return found->second;

	SbColor color(1.0f, 0.0f, 0.0f);
	struct region regionState;
	memset(&regionState, 0, sizeof(regionState));
	regionState.reg_regionid = regionId;
	db_mater_color_region(this->dbip, &regionState);
	if (regionState.reg_mater.ma_color_valid) {
	    color = SbColor(
		static_cast<float>(regionState.reg_mater.ma_color[0]),
		static_cast<float>(regionState.reg_mater.ma_color[1]),
		static_cast<float>(regionState.reg_mater.ma_color[2]));
	}
	this->regionColors.emplace(regionId, color);
	return color;
    }

    struct db_i *dbip;
    std::unordered_map<struct directory *, BObolMaterialObjectState>
	objectStates;
    std::unordered_map<std::string, BObolMaterialPathState> pathStates;
    std::unordered_map<int, SbColor> regionColors;
    BObolMaterialPathNode pathRoot;
};

} // namespace

SbBool
bobol_database_source_fullpath_material_color(
    struct db_i *dbip,
    const struct db_full_path *pathp,
    SbColor &color)
{
    BObolMaterialColorSweep sweep(dbip);
    BObolMaterialPathState resolved;
    if (!sweep.resolve(pathp, resolved))
	return FALSE;
    color = resolved.color;
    return TRUE;
}



int
bobol_database_sources_refresh_material_colors(
    SoBRLDatabaseSource *const *sources,
    size_t sourceCount,
    uint32_t materialRevision,
    struct db_i *dbip)
{
    if (sourceCount && !dbip)
	return -1;
    return SoBRLDatabaseSource::refreshMaterialColors(sources, sourceCount,
	materialRevision, dbip, nullptr, nullptr);
}

static const char *
lookup_name_from_path(const SbString &path)
{
    const char *name = path.getString();
    while (*name == '/')
	name++;
    if (!name[0])
	return NULL;
    return name;
}

static std::string
stable_name_from_path(const char *path, int fallbackToPath)
{
    if (!path)
	return std::string();

    SbString stablePath(path);
    const char *name = lookup_name_from_path(stablePath);
    if (name && name[0])
	return std::string(name);

    return fallbackToPath ? std::string(path) : std::string();
}

static std::string
database_lookup_path_from_source_path(const SbString &path)
{
    const char *name = lookup_name_from_path(path);
    if (!name)
	return std::string();

    std::string lookup;
    for (const char *cp = name; *cp; cp++) {
	if (*cp == '@') {
	    while (cp[1] && cp[1] != '/')
		cp++;
	    continue;
	}
	lookup.push_back(*cp);
    }
    return lookup;
}

const char *
database_source_skip_leading_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

SbString
record_identity_with_revision(const char *identity, uint32_t revision)
{
    SbString ret = identity ? identity : "";
    char revisionString[64] = {0};
    snprintf(revisionString, sizeof(revisionString), "#%u", revision);
    ret += revisionString;
    return ret;
}



int
source_record_draw_mode(const SoBRLDatabaseSource *source)
{
    return source ? source->getEffectiveLodDrawMode() :
	BOBOL_LOD_DRAW_WIRE;
}

SbBool
SoBRLDatabaseSource::usesMeshRealization(void) const
{
    const int roleFlags = this->realizationRoleFlags.getValue();
    if (roleFlags & SoBRLDatabaseSource::REALIZATION_ROLE_MESH)
	return TRUE;

    /*
     * The role flag is an orchestration hint, not the geometry contract.
     * A cached structural publication may temporarily mark the source
     * EXTERNAL before the deferred worker is cloned.  On a fast warm draw
     * that used to make the worker choose legacy wire realization even
     * though the copied view policy explicitly enabled mesh LoD, allowing its
     * final result to overwrite request-bearing PoP entries.
     *
     * Wire BoTs under an active mesh-LoD policy are intrinsically mesh
     * realization: wire edges are derived from the selected triangle prefix.
     */
    if (this->drawMode.getValue() == SoBRLDatabaseSource::WIREFRAME &&
	this->realizationViewDependent.getValue() &&
	this->realizationMeshLodEnabled.getValue() &&
	this->lodBotThreshold.getValue() > 0)
	return TRUE;

    const int representation = this->representationMode.getValue();
    if (representation == SoBRLDatabaseSource::REPRESENTATION_HIDDEN_LINE ||
	representation == SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS)
	return TRUE;

    return this->drawMode.getValue() == SoBRLDatabaseSource::SHADED ?
	TRUE : FALSE;
}

static int
source_uses_evaluated_wire_realization(const SoBRLDatabaseSource *source)
{
    return source && source->representationMode.getValue() ==
	   SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE;
}

static int
source_uses_evaluated_points_realization(const SoBRLDatabaseSource *source)
{
    return source && source->representationMode.getValue() ==
	   SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS;
}

static int
source_uses_evaluated_path_realization(const SoBRLDatabaseSource *source)
{
    return source_uses_evaluated_wire_realization(source) ||
	   source_uses_evaluated_points_realization(source);
}

SbString
source_effective_instance_key(const SoBRLDatabaseSource *source)
{
    if (!source)
	return "";
    const SbString key = source->instanceKey.getValue();
    if (key.getLength() > 0)
	return key;
    return source->path.getValue();
}

static uint64_t
source_stable_compact_handle_id(const SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;
    const uintptr_t databaseIdentity = reinterpret_cast<uintptr_t>(
	source->getDatabase());
    const SbString key = source_effective_instance_key(source);
    const std::pair<uintptr_t, std::string> exactKey(
	databaseIdentity, key.getString());

    /* Public compact handles need a stable scalar, but a digest cannot be
     * their authorization boundary.  Intern the exact database/path tuple
     * and allocate a process-lifetime token which is never recycled. */
    static std::mutex registryMutex;
    static std::map<std::pair<uintptr_t, std::string>, uint64_t> registry;
    static uint64_t nextIdentity = 1;
    std::lock_guard<std::mutex> lock(registryMutex);
    const auto found = registry.find(exactKey);
    if (found != registry.end())
	return found->second;
    const uint64_t identity = bobol_nonzero_identity_take(nextIdentity);
    registry.emplace(exactKey, identity);
    return identity;
}

SbString
source_effective_representation_key(const SoBRLDatabaseSource *source)
{
    if (!source)
	return "";
    const SbString key = source->representationKey.getValue();
    if (key.getLength() > 0)
	return key;
    return source_effective_instance_key(source);
}

static int
source_instance_matches_record_path(const SbString &instanceKey,
				    const char *recordPath)
{
    if (!recordPath || !recordPath[0])
	return 0;
    if (instanceKey.getLength() == 0)
	return 0;

    const char *key = instanceKey.getString();
    if (bu_strcmp(key, recordPath) == 0)
	return 1;
    return bu_strcmp(key[0] == '/' ? key + 1 : key,
		  recordPath[0] == '/' ? recordPath + 1 : recordPath) == 0;
}

static SbString
source_record_identity(const SoBRLDatabaseSource *source,
		       const char *recordPath)
{
    const char *path = recordPath ? recordPath : "";
    if (!source)
	return path;

    const SbString instanceKey = source_effective_instance_key(source);
    const SbString representationKey =
	source_effective_representation_key(source);
    if (source_instance_matches_record_path(instanceKey, path)) {
	if (representationKey.getLength() == 0 ||
	    source_instance_matches_record_path(representationKey, path))
	    return path;
	SbString identity = representationKey;
	if (path[0]) {
	    identity += "::";
	    identity += path;
	}
	return identity;
    }
    if (instanceKey.getLength() == 0 && representationKey.getLength() == 0)
	return path;

    SbString identity = representationKey.getLength() > 0 ?
			representationKey : instanceKey;
    if (path[0]) {
	identity += "::";
	identity += path;
    }
    return identity;
}

static SbString
source_realization_identity(const SoBRLDatabaseSource *source, uint32_t revision)
{
    if (!source)
	return "";

    SbString identity;
    identity.sprintf("dbsource:%s@%s#repr=%s;repr_mode=%d;mode=%d;source=%u;inputs=%u;view=%u;view_dep=%d;csg_lod=%d;mesh_lod=%d;view_scale=%.9g;lod_scale=%.9g;view_dims=%dx%d;lod=%u;curve=%.9g;point=%.9g;abs=%.9g;rel=%.9g;norm=%.9g",
		     source_effective_instance_key(source).getString(),
		     source->path.getValue().getString(),
		     source_effective_representation_key(source).getString(),
		     source->representationMode.getValue(),
		     source->drawMode.getValue(),
		     revision,
		     source->inputsRevision.getValue(),
		     source->viewRevision.getValue(),
		     source->realizationViewDependent.getValue() ? 1 : 0,
		     source->realizationCsgLodEnabled.getValue() ? 1 : 0,
		     source->realizationMeshLodEnabled.getValue() ? 1 : 0,
		     static_cast<double>(source->realizationViewScale.getValue()),
		     static_cast<double>(source->realizationLodScale.getValue()),
		     source->realizationViewWidth.getValue(),
		     source->realizationViewHeight.getValue(),
		     source->lodBotThreshold.getValue(),
		     static_cast<double>(source->realizationCurveScale.getValue()),
		     static_cast<double>(source->realizationPointScale.getValue()),
		     static_cast<double>(source->tessellationAbsTol.getValue()),
		     static_cast<double>(source->tessellationRelTol.getValue()),
		     static_cast<double>(source->tessellationNormTol.getValue()));
    return identity;
}

static SbString
source_realization_identity(const SoBRLDatabaseSource *source)
{
    return source_realization_identity(source, source ? source->sourceRevision.getValue() : 0);
}

static void
convert_vlist(std::vector<SbVec3f> &points, std::vector<int32_t> &commands, const struct bu_list *vhead)
{
    rt_vlist *vp = NULL;

    BU_LIST_EACH(vhead, vp, rt_vlist) {
	for (size_t i = 0; i < vp->nused; i++) {
	    int cmd = -1;
	    switch (vp->cmd[i]) {
		case RT_VLIST_LINE_MOVE:
		case RT_VLIST_POLY_MOVE:
		case RT_VLIST_TRI_MOVE:
		    cmd = SoBRLVListShape::MOVE;
		    break;
		case RT_VLIST_LINE_DRAW:
		case RT_VLIST_POLY_DRAW:
		case RT_VLIST_POLY_END:
		case RT_VLIST_TRI_DRAW:
		case RT_VLIST_TRI_END:
		    cmd = SoBRLVListShape::DRAW;
		    break;
		case RT_VLIST_POINT_DRAW:
		    cmd = SoBRLVListShape::POINT;
		    break;
		default:
		    break;
	    }
	    if (cmd >= 0) {
		points.push_back(SbVec3f(static_cast<float>(vp->pt[i][0]),
					 static_cast<float>(vp->pt[i][1]),
					 static_cast<float>(vp->pt[i][2])));
		commands.push_back(cmd);
	    }
	}
    }
}

/* rt_obj_plot allocates from the process-global rt_vlfree list, whose macros
 * are explicitly non-parallel.  Detached realization and a foreground draw
 * may otherwise corrupt that list between plot, conversion, and return. */
static std::mutex database_source_rt_vlist_mutex;

static int
plot_internal_to_vlist_geometry(
	std::vector<SbVec3f> &points,
	std::vector<int32_t> &commands,
	struct rt_db_internal *intern,
	const struct bg_tess_tol *ttol,
	const struct bn_tol *tol)
{
    if (!intern)
	return -1;

    std::lock_guard<std::mutex> guard(database_source_rt_vlist_mutex);
    struct bu_list vhead;
    BU_LIST_INIT(&vhead);
    const auto release = [](bu_list *head) { RT_FREE_VLIST(&rt_vlfree, head); };
    std::unique_ptr<bu_list, decltype(release)> plotted(&vhead, release);
    int ret = 0;
    {
	BObolPerformanceTimer timer(BOBOL_PERF_PLOT_US);
	if (timer.active())
	    bobol_performance_counter_add(BOBOL_PERF_PLOT_CALLS, 1);
	ret = rt_obj_plot(&vhead, intern, ttol, tol);
    }
    if (ret >= 0) {
	BObolPerformanceTimer timer(BOBOL_PERF_VLIST_CONVERT_US);
	if (timer.active())
	    bobol_performance_counter_add(BOBOL_PERF_VLIST_CONVERT_CALLS, 1);
	convert_vlist(points, commands, &vhead);
	if (intern->idb_type == ID_ANNOT) {
	    const auto *annotation = static_cast<const rt_annot_internal *>(intern->idb_ptr);
	    RT_ANNOT_CK_MAGIC(annotation);
	    /* Model-space plots already include their anchor and plane basis.
	     * The legacy screen-space fallback needs its display anchor restored
	     * after conversion discards the display-matrix commands. */
	    if (!(annotation->flags & RT_ANNOT_MODEL_SPACE)) {
		const SbVec3f anchor(float(annotation->V[X]), float(annotation->V[Y]),
		    float(annotation->V[Z]));
		for (SbVec3f &point : points)
		    point += anchor;
	    }
	}
	if (!points.empty())
	    bobol_performance_counter_add(BOBOL_PERF_VLIST_POINTS,
		static_cast<uint64_t>(points.size()));
    }
    return ret;
}

static SoBRLVListShape *
vlist_from_bot_wireframe(const struct rt_bot_internal *bot)
{
    if (!bot || !bot->vertices || !bot->faces ||
	bot->num_vertices == 0 || bot->num_faces == 0 ||
	bot->num_faces > static_cast<size_t>(INT_MAX / 4))
	return NULL;
    RT_BOT_CK_MAGIC(bot);

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    points.reserve(bot->num_faces * 4);
    commands.reserve(bot->num_faces * 4);

    for (size_t i = 0; i < bot->num_faces; i++) {
	const int *face = &bot->faces[i * 3];
	if (face[0] < 0 || face[1] < 0 || face[2] < 0 ||
	    static_cast<size_t>(face[0]) >= bot->num_vertices ||
	    static_cast<size_t>(face[1]) >= bot->num_vertices ||
	    static_cast<size_t>(face[2]) >= bot->num_vertices)
	    continue;

	const fastf_t *p0 = &bot->vertices[face[0] * 3];
	const fastf_t *p1 = &bot->vertices[face[1] * 3];
	const fastf_t *p2 = &bot->vertices[face[2] * 3];

	points.push_back(SbVec3f(static_cast<float>(p0[X]),
				 static_cast<float>(p0[Y]),
				 static_cast<float>(p0[Z])));
	commands.push_back(SoBRLVListShape::MOVE);
	points.push_back(SbVec3f(static_cast<float>(p1[X]),
				 static_cast<float>(p1[Y]),
				 static_cast<float>(p1[Z])));
	commands.push_back(SoBRLVListShape::DRAW);
	points.push_back(SbVec3f(static_cast<float>(p2[X]),
				 static_cast<float>(p2[Y]),
				 static_cast<float>(p2[Z])));
	commands.push_back(SoBRLVListShape::DRAW);
	points.push_back(SbVec3f(static_cast<float>(p0[X]),
				 static_cast<float>(p0[Y]),
				 static_cast<float>(p0[Z])));
	commands.push_back(SoBRLVListShape::DRAW);
    }

    if (points.empty() || points.size() != commands.size())
	return NULL;

    bobol_performance_counter_add(BOBOL_PERF_VLIST_POINTS,
	static_cast<uint64_t>(points.size()));

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->setLineSet(points.data(), commands.data(),
		      static_cast<int>(points.size()));
    return shape;
}

static int
cad_wire_part_geometry_from_line_set(const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &commands, Obol::PartGeometryBuilder &geometry)
{
    const size_t count = std::min(points.size(), commands.size());
    if (!count)
	return 0;

    Obol::WireRep wire;
    wire.bounds.makeEmpty();
    wire.segmentPoints.reserve(count * 2u);
    wire.segmentIds.reserve(count);
    Obol::PointRep pointRep;
    pointRep.bounds.makeEmpty();
    pointRep.positions.reserve(count);
    pointRep.pointIds.reserve(count);
    bool haveLast = false;
    size_t lastIndex = 0;
    uint32_t segmentIndex = 0;
    for (size_t i = 0; i < count; i++) {
	if (!std::isfinite(points[i][0]) || !std::isfinite(points[i][1]) ||
	    !std::isfinite(points[i][2])) {
	    if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
		bu_log("[obol-timing] rejected non-finite wire point %zu: "
		       "%.9g %.9g %.9g\n", i,
		       static_cast<double>(points[i][0]),
		       static_cast<double>(points[i][1]),
		       static_cast<double>(points[i][2]));
	    return 0;
	}
	const int command = commands[i];
	if (command == SoBRLVListShape::POINT) {
	    const SbVec3f &point = points[i];
	    pointRep.positions.push_back(point);
	    pointRep.pointIds.push_back(static_cast<uint32_t>(i));
	    pointRep.colorValid.push_back(0u);
	    pointRep.colors.push_back(SbColor(1.0f, 1.0f, 1.0f));
	    pointRep.scaleValid.push_back(0u);
	    pointRep.scales.push_back(0.0f);
	    pointRep.normalValid.push_back(0u);
	    pointRep.normals.push_back(SbVec3f(0.0f, 0.0f, 1.0f));
	    pointRep.bounds.extendBy(point);
	    continue;
	}
	if (command == SoBRLVListShape::MOVE) {
	    lastIndex = i;
	    haveLast = true;
	    continue;
	}
	if (command != SoBRLVListShape::DRAW || !haveLast)
	    continue;
	const SbVec3f &a = points[lastIndex];
	const SbVec3f &b = points[i];
	wire.segmentPoints.push_back(a);
	wire.segmentPoints.push_back(b);
	wire.segmentIds.push_back(segmentIndex++);
	wire.bounds.extendBy(a);
	wire.bounds.extendBy(b);
	lastIndex = i;
    }
    if (wire.segmentPoints.empty() && pointRep.positions.empty())
	return 0;
    if (!wire.segmentPoints.empty())
	geometry.wire = std::move(wire);
    if (!pointRep.positions.empty())
	geometry.points = std::move(pointRep);
    return 1;
}

static int
cad_wire_part_geometry_from_aabb(const SbBox3f &bounds,
	Obol::PartGeometryBuilder &geometry)
{
    if (bounds.isEmpty())
	return 0;

    const SbVec3f bmin = bounds.getMin();
    const SbVec3f bmax = bounds.getMax();
    const SbVec3f corners[8] = {
	SbVec3f(bmin[0], bmin[1], bmin[2]),
	SbVec3f(bmax[0], bmin[1], bmin[2]),
	SbVec3f(bmax[0], bmax[1], bmin[2]),
	SbVec3f(bmin[0], bmax[1], bmin[2]),
	SbVec3f(bmin[0], bmin[1], bmax[2]),
	SbVec3f(bmax[0], bmin[1], bmax[2]),
	SbVec3f(bmax[0], bmax[1], bmax[2]),
	SbVec3f(bmin[0], bmax[1], bmax[2])
    };
    static const int edges[12][2] = {
	{0, 1}, {1, 2}, {2, 3}, {3, 0},
	{4, 5}, {5, 6}, {6, 7}, {7, 4},
	{0, 4}, {1, 5}, {2, 6}, {3, 7}
    };

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    points.reserve(24);
    commands.reserve(24);
    for (size_t i = 0; i < 12; i++) {
	points.push_back(corners[edges[i][0]]);
	commands.push_back(SoBRLVListShape::MOVE);
	points.push_back(corners[edges[i][1]]);
	commands.push_back(SoBRLVListShape::DRAW);
    }
    if (!cad_wire_part_geometry_from_line_set(points, commands, geometry))
	return 0;
    /* This helper is used only for the mesh LoD AABB threshold path. */
    geometry.subpixelProxyEligible = true;
    geometry.structuralProxy = true;
    return 1;
}

static int
cad_source_mesh_request_from_bot(BObolSourceMeshRequest &request,
	const struct rt_bot_internal *bot)
{
    request.clear();
    if (!bot || !bot->vertices || !bot->faces || !bot->num_vertices ||
	!bot->num_faces)
	return 0;
    RT_BOT_CK_MAGIC(bot);

    SbBox3f bounds;
    bounds.makeEmpty();
    for (size_t i = 0; i < bot->num_vertices; i++) {
	const fastf_t *vertex = &bot->vertices[i * 3u];
	if (!isfinite(vertex[X]) || !isfinite(vertex[Y]) ||
	    !isfinite(vertex[Z]))
	    return 0;
	bounds.extendBy(SbVec3f(static_cast<float>(vertex[X]),
	    static_cast<float>(vertex[Y]), static_cast<float>(vertex[Z])));
    }
    if (bounds.isEmpty())
	return 0;
    request.faceCount = static_cast<uint64_t>(bot->num_faces);
    request.pointCount = static_cast<uint64_t>(bot->num_vertices);
    request.bounds = bounds;
    request.meshAssetBounds = bounds;
    return 1;
}


static double
nonnegative_or_default(float value, double defaultValue)
{
    return value >= 0.0f ? static_cast<double>(value) : defaultValue;
}

static struct bg_tess_tol
source_tess_tol(const SoBRLDatabaseSource *source)
{
    struct bg_tess_tol ttol = BG_TESS_TOL_INIT_TOL;

    if (source) {
	ttol.abs = nonnegative_or_default(source->tessellationAbsTol.getValue(), ttol.abs);
	ttol.rel = nonnegative_or_default(source->tessellationRelTol.getValue(), ttol.rel);
	ttol.norm = nonnegative_or_default(source->tessellationNormTol.getValue(), ttol.norm);
    }

    return ttol;
}

static SbBool
source_view_lod_active(const SoBRLDatabaseSource *source)
{
    return source && source->realizationViewDependent.getValue() ? TRUE : FALSE;
}

static SbBool
source_csg_lod_active(const SoBRLDatabaseSource *source)
{
    return source_view_lod_active(source) &&
	   source->realizationCsgLodEnabled.getValue() ? TRUE : FALSE;
}

static SbBool
source_mesh_lod_active(const SoBRLDatabaseSource *source)
{
    return source_view_lod_active(source) &&
	   source->realizationMeshLodEnabled.getValue() ? TRUE : FALSE;
}

static void
source_view_info(struct bv_view_info *info,
		 const SoBRLDatabaseSource *source)
{
    bv_view_info_init(info);
    if (!info || !source)
	return;

    info->size = source->realizationViewScale.getValue();
    info->width = source->realizationViewWidth.getValue();
    info->height = source->realizationViewHeight.getValue();
    info->lod.scale = source->realizationLodScale.getValue();
    info->lod.curve_scale = source->realizationCurveScale.getValue();
    info->lod.point_scale = source->realizationPointScale.getValue();
    info->lod.bot_threshold = source->realizationBotThreshold.getValue();
    bv_view_info_sanitize(info);
}

static fastf_t
source_lod_solid_size(const SoBRLDatabaseSource *source,
		      const SbBox3f &localBounds)
{
    if (source && source->drawSizeValid.getValue() &&
	source->drawSize.getValue() > 0.0f)
	return source->drawSize.getValue();

    if (!localBounds.isEmpty()) {
	const SbVec3f bmin = localBounds.getMin();
	const SbVec3f bmax = localBounds.getMax();
	return (bmax - bmin).length();
    }

    if (source && source->realizationViewScale.getValue() > 0.0f)
	return source->realizationViewScale.getValue();

    return 1.0;
}

static void
source_lod_cache_key_append(std::string &cacheKey,
			    const SoBRLDatabaseSource *source,
			    const SbBox3f &localBounds,
			    bool viewIndependent = false)
{
    /* Normal wire geometry is source-local and independent of viewport
     * policy.  Only CSG LoD realizes a view-dependent line payload here;
     * mesh/BOT proxy variants carry their geometry-affecting threshold in
     * their own key. */
    if (viewIndependent) {
	/* Canonical BREP wire and shaded assets never encode a camera.  They do
	 * encode the source-space tessellation contract, so changing an adaptive
	 * tessellation setting cannot accidentally reuse an in-memory payload
	 * produced for a different error bound. */
	char suffix[256] = {0};
	snprintf(suffix, sizeof(suffix),
	    ":brep-source-lod:%s:%d:%.17g:%.17g:%.17g",
	    BOBOL_MESH_LOD_PROVIDER_VERSION,
	    source ? source->drawMode.getValue() :
		SoBRLDatabaseSource::WIREFRAME,
	    static_cast<double>(source ?
		source->tessellationAbsTol.getValue() : -1.0f),
	    static_cast<double>(source ?
		source->tessellationRelTol.getValue() : -1.0f),
	    static_cast<double>(source ?
		source->tessellationNormTol.getValue() : -1.0f));
	cacheKey += suffix;
	return;
    }
    if (!source_csg_lod_active(source))
	return;

    char suffix[256] = {0};
    snprintf(suffix, sizeof(suffix),
	     ":view-lod:%d:%d:%u:%.9g:%.9g:%dx%d:%.9g:%.9g:%.9g",
	     source->realizationCsgLodEnabled.getValue() ? 1 : 0,
	     source->realizationMeshLodEnabled.getValue() ? 1 : 0,
	     source->realizationBotThreshold.getValue(),
	     static_cast<double>(source->realizationViewScale.getValue()),
	     static_cast<double>(source->realizationLodScale.getValue()),
	     source->realizationViewWidth.getValue(),
	     source->realizationViewHeight.getValue(),
	     static_cast<double>(source->realizationCurveScale.getValue()),
	     static_cast<double>(source->realizationPointScale.getValue()),
	     static_cast<double>(source_lod_solid_size(source, localBounds)));
    cacheKey += suffix;
}

static void
primitive_realization_line_set_free(
    struct rt_primitive_lod_realization *realization);

static int32_t
primitive_realization_command_to_vlist_command(int command);

static SoBRLVListShape *
vlist_from_primitive_realization_line_set(
    struct rt_primitive_lod_realization *realization,
    const char *geometryKind);

static SbBool
local_bounds_from_internal(struct rt_db_internal *intern,
			   SbBox3f &bounds);

static SbBool
node_is_auxiliary_vlist(const SoNode *node)
{
    if (!node || !node->isOfType(SoBRLVListShape::getClassTypeId()))
	return FALSE;

    const SoBRLVListShape *shape = static_cast<const SoBRLVListShape *>(node);
    return bu_strcmp(shape->recordRole.getValue().getString(), "auxiliary") == 0 ?
	   TRUE : FALSE;
}

static SbBool
node_is_auxiliary_source(const SoNode *node)
{
    if (!node || !node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	return FALSE;

    const SoBRLDatabaseSource *source =
	static_cast<const SoBRLDatabaseSource *>(node);
    return source->auxiliarySource.getValue();
}

static const char *
database_source_placement_transform_name(void)
{
    return "__bobol_source_placement";
}

SbBool
node_is_source_placement_transform(const SoNode *node)
{
    return node &&
	   node->isOfType(SoMatrixTransform::getClassTypeId()) &&
	   node->getName() == SbName(database_source_placement_transform_name()) ?
	   TRUE : FALSE;
}

static SoMatrixTransform *
source_placement_transform(SoBRLDatabaseSource *source)
{
    if (!source)
	return NULL;

    for (int i = 0; i < source->getNumChildren(); i++) {
	SoNode *child = source->getChild(i);
	if (node_is_source_placement_transform(child))
	    return static_cast<SoMatrixTransform *>(child);
    }

    return NULL;
}

static const SoMatrixTransform *
source_placement_transform(const SoBRLDatabaseSource *source)
{
    return source_placement_transform(const_cast<SoBRLDatabaseSource *>(source));
}

static SoMatrixTransform *
ensure_source_placement_transform(SoBRLDatabaseSource *source)
{
    if (!source)
	return NULL;

    SoMatrixTransform *transform = source_placement_transform(source);
    if (!transform) {
	transform = new SoMatrixTransform;
	transform->setName(SbName(database_source_placement_transform_name()));
	source->insertChild(transform, 0);
	return transform;
    }

    int transformIndex = -1;
    for (int i = 0; i < source->getNumChildren(); i++) {
	if (source->getChild(i) == transform) {
	    transformIndex = i;
	    break;
	}
    }
    if (transformIndex > 0) {
	transform->ref();
	source->removeChild(transformIndex);
	source->insertChild(transform, 0);
	transform->unref();
    }

    return transform;
}

static int
remove_source_placement_transform(SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;

    int removed = 0;
    for (int i = source->getNumChildren() - 1; i >= 0; i--) {
	if (node_is_source_placement_transform(source->getChild(i))) {
	    source->removeChild(i);
	    removed = 1;
	}
    }
    return removed;
}

static int
sync_source_placement_transform(SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;

    if (!source->drawMatrixValid.getValue()) {
	return remove_source_placement_transform(source);
    }

    const int hadTransform = source_placement_transform(source) ? 1 : 0;
    SoMatrixTransform *transform = ensure_source_placement_transform(source);
    if (!transform)
	return 0;

    int changed = hadTransform ? 0 : 1;
    const SbMatrix matrix = source->drawMatrix.getValue();
    if (!transform->matrix.getValue().equals(matrix, 0.000001f)) {
	transform->matrix = matrix;
	changed = 1;
    }
    return changed;
}

static void
database_source_add_realized_child(SoBRLDatabaseSource *source,
				   SoNode *child)
{
    if (!source || !child)
	return;

    (void)sync_source_placement_transform(source);
    source->addChild(child);
}

static bool
source_retains_realized_child(SoNode *child)
{
    return node_is_source_placement_transform(child) ||
	node_is_auxiliary_vlist(child) || node_is_auxiliary_source(child);
}

static std::unique_ptr<SoChildList::Removal>
prepare_non_auxiliary_child_removal(SoGroup *group)
{
    if (!group || !group->getNumChildren())
	return nullptr;
    std::vector<int> indices;
    for (int i = 0; i < group->getNumChildren(); ++i)
	if (!source_retains_realized_child(group->getChild(i)))
	    indices.push_back(i);
    return group->getChildren()->prepareRemoval(std::move(indices));
}

static void
remove_non_auxiliary_children(SoGroup *group)
{
    if (!group)
	return;

    for (int i = group->getNumChildren() - 1; i >= 0; i--) {
	SoNode *child = group->getChild(i);
	if (source_retains_realized_child(child))
	    continue;
	group->removeChild(i);
    }
}

/* Line PartGeometry is representation-compatible across the legacy wire walk
 * and the mesh-role wire path used by LoD-enabled views.  Query both line
 * namespaces only from a wire realization; shaded lookups must never consume
 * these payloads. */
static const BObolCachedPartGeometry *
find_wire_cad_geometry_any(
    const BObolDatabaseSourceRealizationCache *cache,
    const std::string &key)
{
    if (!cache)
	return NULL;
    const BObolCachedPartGeometry *geometry =
	cache->findWireCadGeometry(key);
    return geometry ? geometry : cache->findMeshVListCadGeometry(key);
}

static void
cache_mesh_cad_source_request(
    BObolDatabaseSourceRealizationCache *cache,
    const std::string &key, const BObolSourceMeshRequest &request)
{
    if (!cache || key.empty() || !request.faceCount || !request.pointCount ||
	request.bounds.isEmpty())
	return;

    BObolRealizationCacheMap<BObolCachedPartGeometry>::iterator found =
	cache->sharedMeshCadGeometry.find(key);
    if (found == cache->sharedMeshCadGeometry.end())
	return;
    found->second.sourceMeshRequest = request;
    found->second.sourceMeshRequestValid = true;
}

template <typename ShapeT>
static int
cacheable_shared_geometry_key(ShapeT *shape,
			      ShapeT *geometry,
			      std::string &key)
{
    if (!shape || !geometry || geometry == shape)
	return 0;
    if (!geometry->sharedSource.getValue())
	return 0;

    const char *name = geometry->geometryName.getValue().getString();
    if (!name || !name[0])
	return 0;

    key = name;
    return 1;
}

static void
seed_realization_cache_from_node(SoNode *node,
				 BObolDatabaseSourceRealizationCache *cache,
				 int meshRealization)
{
    if (!node || !cache)
	return;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	SoBRLVListShape *geometry = shape->getSharedGeometrySource();
	std::string key;
	if (cacheable_shared_geometry_key(shape, geometry, key)) {
	    if (meshRealization)
		bobol_cache_geometry_reference(cache->sharedMeshVListGeometry, key, geometry);
	    else
		bobol_cache_geometry_reference(cache->sharedWireGeometry, key, geometry);
	}
	return;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	SoBRLMeshShape *geometry = shape->getSharedGeometrySource();
	std::string key;
	if (meshRealization && cacheable_shared_geometry_key(shape,
		geometry, key))
	    bobol_cache_geometry_reference(cache->sharedMeshGeometry, key, geometry);
	return;
    }

    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    seed_realization_cache_from_node(group->getChild(i), cache,
					     meshRealization);
    }
}

static SbBool source_bounds_for_realized_node(const SoNode *node,
	const SbMatrix &matrix, SbBox3f &bounds);

void
bobol_database_source_seed_realization_cache(
    SoBRLDatabaseSource *source,
    BObolDatabaseSourceRealizationCache *cache)
{
    if (!source || !cache ||
	source->realizationStatus.getValue() !=
	SoBRLDatabaseSource::REALIZED ||
	source->needsRealization())
	return;

    const int roleFlags = source->realizationRoleFlags.getValue();
    if (roleFlags & SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL)
	return;

    if (source->hasCompactInstanceIndex()) {
	source->seedCompactRealizationCache(cache);
	return;
    }

    seed_realization_cache_from_node(source, cache,
				     source->usesMeshRealization());
}

void
SoBRLDatabaseSource::seedCompactRealizationCache(
    BObolDatabaseSourceRealizationCache *cache) const
{
    if (!cache || !this->d->compactIndex)
	return;

    const bool meshRealization = this->usesMeshRealization() ? true : false;

    for (const BObolCompactInstanceEntry &entry :
	 this->d->compactIndex->entries) {
	/* This entry was produced for the source's realized camera policy.  The
	 * detached successor builds cache keys from its new policy, so inserting
	 * the old payload under that key would falsely certify stale CSG geometry
	 * as current.  Immutable mesh/PoP and ordinary plotted geometry remain
	 * safe to seed and reuse. */
	if (entry.viewDependentCsgGeometry)
	    continue;
	const char *name = entry.semantic.sourceName.getString();
	if (!name || !name[0])
	    continue;

	const SbBox3f localBounds = compact_part_geometry_bounds(entry.geometry);
	std::string cacheKey(name);
	/* Realization lookup keys are formed before an internal is fetched, so
	 * their LoD size fallback intentionally starts with empty bounds. */
	SbBox3f lookupBounds;
	const bool viewIndependent = BU_STR_EQUAL(
	    entry.shapeSummary.sourceType.getString(), "brep");
	source_lod_cache_key_append(cacheKey, this, lookupBounds,
	    viewIndependent);

	if (entry.geometry) {
	    BObolCachedPartGeometry cached;
	    cached.geometry = entry.geometry;
	    cached.sourceType = entry.shapeSummary.sourceType.getString();
	    cached.geometryKind = entry.shapeSummary.geometryKind.getString();
	    cached.geometryTransform = entry.geometryTransform;
	    cached.bounds = database_source_transform_bounds(localBounds,
		entry.geometryTransform);
	    cached.lodBacked = entry.lodBacked;
	    cached.sourceMeshRequestValid = entry.sourceMeshRequestValid;
	    if (cached.sourceMeshRequestValid)
		cached.sourceMeshRequest = entry.sourceMeshRequest;
	    if (entry.meshGeometry)
		cache->sharedMeshCadGeometry[cacheKey] = cached;
	    else if ((entry.wireGeometry || entry.pointGeometry) && meshRealization)
		cache->sharedMeshVListCadGeometry[cacheKey] = cached;
	    else if (entry.wireGeometry || entry.pointGeometry)
		cache->sharedWireCadGeometry[cacheKey] = cached;
	}
	if (!localBounds.isEmpty() &&
	    (entry.wireGeometry || entry.pointGeometry) && !meshRealization)
	    cache->storeWireBounds(cacheKey, localBounds);
    }
}

struct compact_stream_lod_reuse_entry {
    struct bg_trimesh_pca_signature signature;
    std::array<int64_t, 3> bucket = {};
    fastf_t comparisonTolerance = VUNITIZE_TOL;
    bool signatureValid = false;
    uint64_t sampleFingerprint = 0;
    bool sampleFingerprintValid = false;
    struct directory *dp = NULL;
    std::string cacheKey;
    BObolSourceMeshRequest sourceMeshRequest;
    size_t vertexCount = 0;
    size_t faceCount = 0;
    unsigned char mode = 0;
    unsigned char orientation = 0;
};

struct realize_walk_data {
    realize_walk_data(void) :
	source(NULL),
	cache(NULL),
	revision(0),
	visited_leaves(0),
	realized_shapes(0),
	failed_shapes(0),
	compact_ordinal(0),
	compact_unsupported(0),
	compact_bounds_valid(FALSE),
	stream_sink(NULL),
	material_sweep(NULL),
	stream_lod_cached_representative_dp(NULL),
	stream_lod_cached_representative_valid(false)
    {
	compact_bounds.makeEmpty();
	RT_DB_INTERNAL_INIT(&stream_lod_cached_representative);
    }

    ~realize_walk_data(void)
    {
	if (stream_lod_cached_representative_valid)
	    rt_db_free_internal(&stream_lod_cached_representative);
    }

    SoBRLDatabaseSource *source;
    BObolDatabaseSourceRealizationCache *cache;
    uint32_t revision;
    size_t visited_leaves;
    int realized_shapes;
    int failed_shapes;
    SbString diagnostic;
    std::unique_ptr<BObolCompactInstanceIndex> compact_index;
    /* Optional worker->pump hand-off: when set, each completed compact
     * occurrence is pushed here as it is realized so the progressive pump can
     * stream geometry onto the standing compact root incrementally. */
    int compact_ordinal;
    int compact_unsupported;
    SbBox3f compact_bounds;
    SbBool compact_bounds_valid;
    BObolCompactOccurrenceStream *stream_sink;
    /* One cached path/material resolver per database walk.  Region-table
     * fallback must not re-import every combination for every occurrence. */
    void *material_sweep;
    std::unordered_set<std::string> compact_seen_instances;
    std::unordered_map<std::string, uint32_t> compact_occurrence_counts;
    /* Online transformed-copy representatives for progressive streaming.
     * Keeping signatures and directory identities (not mesh arrays) preserves
     * first-leaf latency while bounding matching memory to one candidate plus
     * one temporarily re-imported representative. */
    std::vector<compact_stream_lod_reuse_entry> stream_lod_reuse;
    size_t stream_lod_asset_hits = 0;
    size_t stream_lod_pca_deferred = 0;
    size_t stream_lod_pca_evaluated = 0;
    size_t stream_lod_pca_reused = 0;
    /* One representative internal is sufficient for the normal transformed-
     * copy run.  Keeping that one import across adjacent candidates has the
     * same peak memory as the old candidate+temporary-representative pair,
     * but avoids rereading a multi-hundred-megabyte BoT for every duplicate.
     */
    struct directory *stream_lod_cached_representative_dp;
    struct rt_db_internal stream_lod_cached_representative;
    bool stream_lod_cached_representative_valid;
    size_t stream_lod_representative_imports = 0;
};

static void
realize_walk_extend_bounds(realize_walk_data *data, const SbBox3f &bounds)
{
    if (!data)
	return;

    if (bounds.isEmpty())
	return;

    data->compact_bounds.extendBy(bounds);
    data->compact_bounds_valid = TRUE;
}

static int cad_vlist_part_geometry(const SoBRLVListShape *shape,
	Obol::PartGeometryBuilder &geometry);
static int cad_mesh_part_geometry(const SoBRLMeshShape *shape,
	Obol::PartGeometryBuilder &geometry);
struct compact_occurrence_build {
    BObolCompactOccurrence occurrence;
    SoBRLCadAssembly::InstanceSemantic semantic;
    Obol::InstanceStyle normalStyle;
    Obol::InstanceStyle selectedStyle;
    Obol::InstanceStyle highlightedStyle;
    SbBool stylesValid = FALSE;
    SbBool dashed = FALSE;
};
static SoBRLCadAssembly::InstanceSemantic compact_semantic_from_summary(
	const BObolRealizedShapeSummary &summary);
static void compact_add_occurrence(SoBRLDatabaseSource *source,
	BObolCompactInstanceIndex &index,
	const compact_occurrence_build &input, int &ordinal,
	int &unsupported,
	const BObolCompactOccurrenceRegistryState *presentation = nullptr);
static SbBool source_bounds_for_realized_node(const SoNode *node,
	const SbMatrix &matrix, SbBox3f &bounds);
static bool compact_stream_lod_transformed_reuse(
	realize_walk_data *data, struct db_i *dbip, struct directory *dp,
	const char *path, const std::string &cacheKey,
	const struct rt_bot_internal *bot,
	const BObolSourceMeshRequest &sourceMeshRequest);

/* Multiple view controllers may detach and realize the same database root at
 * the same time.  The persistent LoD asset map makes the second realization
 * cheap only after the first has published it; without single-flight
 * coordination both workers can import and PCA-match every large BoT in
 * parallel.  Serialize first-use progressive realization per database/root,
 * while unrelated roots and databases remain independent.  Weak entries keep
 * the registry bounded after a flight finishes. */
static std::shared_ptr<std::mutex>
compact_stream_lod_realization_mutex(const struct db_i *dbip,
	const char *treeName)
{
    struct registry_state {
	std::mutex mutex;
	std::unordered_map<std::string, std::weak_ptr<std::mutex>> entries;
	size_t acquisitions = 0;
    };
    /*
     * The global realization coordinator owns process-lifetime worker
     * threads.  A normal function-local registry destructor can therefore
     * run before the coordinator's destructor and race a final in-flight
     * lookup during shared-library shutdown.  This tiny, bounded registry is
     * deliberately process-lifetime storage; expired weak entries are still
     * scavenged during normal operation.
     */
    static registry_state *registry = new registry_state;

    std::string key = dbip && dbip->dbi_filename ?
	dbip->dbi_filename : "<memory>";
    key += '\n';
    key += treeName ? treeName : "";

    std::lock_guard<std::mutex> guard(registry->mutex);
    if ((++registry->acquisitions & 0xffu) == 0) {
	for (auto it = registry->entries.begin();
		it != registry->entries.end();) {
	    if (it->second.expired())
		it = registry->entries.erase(it);
	    else
		++it;
	}
    }
    const auto found = registry->entries.find(key);
    if (found != registry->entries.end()) {
	std::shared_ptr<std::mutex> mutex = found->second.lock();
	if (mutex)
	    return mutex;
    }
    std::shared_ptr<std::mutex> mutex = std::make_shared<std::mutex>();
    registry->entries[key] = mutex;
    return mutex;
}

/* Hand one just-realized occurrence to the progressive pump for incremental
 * streaming, if a sink is attached and the job has not been cancelled. */
static inline void
realize_walk_stream_push(struct realize_walk_data *data,
	const BObolCompactOccurrence &occurrence)
{
    if (data && data->stream_sink && occurrence.geometry &&
	!data->stream_sink->isCancelled())
	data->stream_sink->push(occurrence);
}

static std::string
realize_walk_instance_identity(const struct db_tree_state *tsp,
	const struct db_full_path *pathp)
{
    std::string key;
    if (!tsp || !pathp)
	return key;
    const uint32_t version = 1;
    const int sofar = tsp->ts_sofar &
	(TS_SOFAR_MINUS | TS_SOFAR_INTER | TS_SOFAR_REGION);
    auto append = [&key](const void *data, size_t size) {
	key.append(static_cast<const char *>(data), size);
    };
    append(&version, sizeof(version));
    append(&pathp->fp_len, sizeof(pathp->fp_len));
    append(&sofar, sizeof(sofar));
    append(tsp->ts_mat, sizeof(tsp->ts_mat));
    for (size_t i = 0; i < pathp->fp_len; i++) {
	const struct directory *dp = pathp->fp_names ? pathp->fp_names[i] : NULL;
	const char *name = dp && dp->d_namep ? dp->d_namep : "";
	const size_t len = strlen(name) + 1;
	const int cinst = pathp->fp_cinst ?
	    DB_FULL_PATH_GET_COMB_INST(pathp, i) : 0;
	append(&len, sizeof(len));
	append(name, len);
	append(&cinst, sizeof(cinst));
    }
    return key;
}

static std::string
realize_walk_occurrence_identity(const struct db_full_path *pathp)
{
    std::string key;
    if (!pathp)
	return key;
    const uint32_t version = 1;
    auto append = [&key](const void *data, size_t size) {
	key.append(static_cast<const char *>(data), size);
    };
    append(&version, sizeof(version));
    append(&pathp->fp_len, sizeof(pathp->fp_len));
    for (size_t i = 0; i < pathp->fp_len; i++) {
	const struct directory *dp = pathp->fp_names ? pathp->fp_names[i] : NULL;
	const char *name = dp && dp->d_namep ? dp->d_namep : "";
	const size_t len = strlen(name) + 1;
	const int cinst = pathp->fp_cinst ?
	    DB_FULL_PATH_GET_COMB_INST(pathp, i) : 0;
	append(&len, sizeof(len));
	append(name, len);
	append(&cinst, sizeof(cinst));
    }
    return key;
}

static void
compact_apply_walk_identity(const SoBRLDatabaseSource *source,
    BObolCompactInstanceIndex &index,
    size_t previousEntryCount, const struct db_tree_state *tsp,
    const struct db_full_path *pathp, const std::string &walkIdentity,
    uint32_t duplicateOrdinal)
{
    if (!tsp || !pathp || index.entries.size() <= previousEntryCount ||
	index.instances.empty())
	return;

    BObolCompactInstanceEntry &entry = index.entries.back();
    std::string occurrenceKey = walkIdentity;
    if (!occurrenceKey.empty()) {
	const SbString sourceKey = source_effective_instance_key(source);
	occurrenceKey.insert(0, sourceKey.getString(), sourceKey.getLength());
	if (duplicateOrdinal > 0) {
	    occurrenceKey.push_back('\0');
	    for (int i = 0; i < 4; i++)
		occurrenceKey.push_back(static_cast<char>(
		    duplicateOrdinal >> (i * 8)));
	}
	const Obol::InstanceId instance =
	    Obol::CadIdBuilder::instanceId(occurrenceKey);
	entry.instance = instance;
	index.instances.back().instance = instance;
	/* db_path_to_string omits a direct combination member's occurrence
	 * number.  Preserve the first legacy path and distinguish subsequent
	 * identical paths so an exact-path edit remains occurrence-local. */
	if (duplicateOrdinal > 0) {
	    SbString occurrencePath;
	    occurrencePath.sprintf("%s@%u", entry.semantic.path.getString(),
		duplicateOrdinal);
	    entry.semantic.path = occurrencePath;
	    entry.shapeSummary.path = occurrencePath;
	    index.instances.back().record.childName = occurrencePath.getString();
	}
	char key[96] = {0};
	snprintf(key, sizeof(key), "compact:%016llx:%016llx",
	    static_cast<unsigned long long>(instance.w0),
	    static_cast<unsigned long long>(instance.w1));
	entry.instanceKey = key;
	entry.semantic.sourceInstanceKey = key;
    }
    entry.occurrenceIndex = pathp->fp_cinst && pathp->fp_len ?
	static_cast<uint32_t>(DB_FULL_PATH_GET_COMB_INST(pathp,
		pathp->fp_len - 1)) : 0;
    entry.booleanOperation = (tsp->ts_sofar & TS_SOFAR_MINUS) ?
	SoBRLDatabaseSource::BOOLEAN_SUBTRACT :
	((tsp->ts_sofar & TS_SOFAR_INTER) ?
	 SoBRLDatabaseSource::BOOLEAN_INTERSECT :
	 SoBRLDatabaseSource::BOOLEAN_UNION);
    Obol::InstanceRecord &record = index.instances.back().record;
    record.occurrenceIndex = entry.occurrenceIndex;
    record.boolOp = entry.booleanOperation ==
	SoBRLDatabaseSource::BOOLEAN_SUBTRACT ? 1 :
	(entry.booleanOperation == SoBRLDatabaseSource::BOOLEAN_INTERSECT ? 2 : 0);
}

/* compact_apply_walk_identity may refine the semantic path (notably for
 * repeated members) and the boolean/occurrence identity after the occurrence
 * has been appended.  Stream that authoritative identity, rather than the
 * pre-append input copy, or the live registry and the completed index can
 * disagree about which CAD occurrence owns a payload. */
static inline void
realize_walk_stream_push_current(struct realize_walk_data *data,
	const BObolCompactOccurrence &input,
	const BObolCompactInstanceIndex &index, size_t previousEntryCount)
{
    if (index.entries.size() <= previousEntryCount)
	return;

    const BObolCompactInstanceEntry &entry = index.entries.back();
    BObolCompactOccurrence occurrence = input;
    occurrence.summary.path = entry.semantic.path;
    occurrence.summary.sourceName = entry.semantic.sourceName;
    occurrence.occurrenceIndex = entry.occurrenceIndex;
    occurrence.booleanOperation = entry.booleanOperation;
    realize_walk_stream_push(data, occurrence);
}

static void
canonicalize_affine_tail(SbMatrix &matrix)
{
    /* BRL-CAD object and draw transforms are affine by contract.  Relative
     * frame calculations may nevertheless leave roundoff in the homogeneous
     * column; preserving that residue turns a valid CAD placement into a
     * projective matrix at the retained-scene validation boundary. */
    matrix[0][3] = 0.0f;
    matrix[1][3] = 0.0f;
    matrix[2][3] = 0.0f;
    matrix[3][3] = 1.0f;
}

static SbMatrix
mat_to_sbmatrix(const mat_t mat)
{

    SbMatrix matrix(
	       static_cast<float>(mat[0]), static_cast<float>(mat[4]),
	       static_cast<float>(mat[8]), static_cast<float>(mat[12]),
	       static_cast<float>(mat[1]), static_cast<float>(mat[5]),
	       static_cast<float>(mat[9]), static_cast<float>(mat[13]),
	       static_cast<float>(mat[2]), static_cast<float>(mat[6]),
	       static_cast<float>(mat[10]), static_cast<float>(mat[14]),
	       static_cast<float>(mat[3]), static_cast<float>(mat[7]),
	       static_cast<float>(mat[11]), static_cast<float>(mat[15]));
    canonicalize_affine_tail(matrix);
    return matrix;
}

static std::string database_source_full_path_string(const SbString &path);

static SoSeparator *
realize_matrix_leaf_separator(const SbMatrix &matrix)
{
    SoSeparator *leaf = new SoSeparator;
    if (matrix.equals(SbMatrix::identity(), 0.0f))
	return leaf;

    SoMatrixTransform *transform = new SoMatrixTransform;
    transform->matrix = matrix;
    leaf->addChild(transform);
    return leaf;
}

static SoSeparator *
realize_instance_leaf_separator(const struct db_tree_state *tsp)
{
    return realize_matrix_leaf_separator(tsp ? mat_to_sbmatrix(tsp->ts_mat) :
	SbMatrix::identity());
}

/* The direct-leaf fast path imports the terminal primitive without walking
 * its parent combinations.  Preserve the path placement explicitly; semantic
 * bounds are already evaluated in model coordinates, and publishing local
 * geometry with an identity placement would put the representation outside
 * the camera selected from those bounds.  If the path cannot be resolved,
 * decline the optimization and let the ordinary tree walk diagnose it. */
static SbBool
database_source_path_matrix(struct db_i *dbip, const SbString &sourcePath,
	SbMatrix &matrix)
{
    matrix.makeIdentity();
    if (!dbip)
	return FALSE;

    const std::string fullPath =
	database_source_full_path_string(sourcePath);
    if (fullPath.empty())
	return FALSE;

    struct db_full_path path;
    db_full_path_init(&path);
    if (db_string_to_path(&path, dbip, fullPath.c_str()) != 0 ||
	path.fp_len == 0 || path.fp_len > static_cast<size_t>(INT_MAX)) {
	db_free_full_path(&path);
	return FALSE;
    }

    mat_t pathMatrix;
    MAT_IDN(pathMatrix);
    const int valid = db_path_to_mat(dbip, &path, pathMatrix,
	static_cast<int>(path.fp_len) - 1);
    db_free_full_path(&path);
    if (!valid)
	return FALSE;

    matrix = mat_to_sbmatrix(pathMatrix);
    return TRUE;
}

static std::string
realize_geometry_cache_key(const struct directory *dp)
{
    return (dp && dp->d_namep) ? std::string(dp->d_namep) : std::string();
}

static void
set_walk_diagnostic(struct realize_walk_data *data,
		    const struct db_full_path *pathp,
		    const char *reason)
{
    if (!data || data->diagnostic.getLength() > 0)
	return;

    char *path = pathp ? db_path_to_string(pathp) : NULL;
    SbString msg;
    msg.sprintf("%s: %s", path ? path : "<unknown>", reason ? reason : "realization failed");
    data->diagnostic = msg;
    if (path)
	bu_free(path, "db_path_to_string");
}

static const char *
primitive_type_label(const struct rt_db_internal *intern)
{
    if (!intern || !intern->idb_meth || !intern->idb_meth->ft_label[0])
	return "unknown";
    return intern->idb_meth->ft_label;
}

static bool
primitive_is_annotation(int internalType, const char *typeLabel)
{
    if (internalType == ID_ANNOT)
	return true;
    return typeLabel && (BU_STR_EQUAL(typeLabel, "annot") ||
			 BU_STR_EQUAL(typeLabel, "annotation"));
}

/* Some database primitives do not have a finite shaded surface of their own.
 * They remain valid members of a shaded display, but must use their plotting
 * representation instead of making the entire aggregate mesh realization
 * fail.  A half-space is the important case: its infinite surface cannot be
 * tessellated as an isolated leaf. */
static bool
primitive_uses_wire_in_mesh_mode(int internalType)
{
    return internalType == ID_HALF || internalType == ID_SKETCH ||
	   internalType == ID_ANNOT;
}

static uint32_t
internal_payload_magic(const struct rt_db_internal *intern)
{
    if (!intern || !intern->idb_ptr)
	return 0;
    return *((const uint32_t *)intern->idb_ptr);
}

static bool
internal_payload_magic_valid(const struct rt_db_internal *intern)
{
    if (!intern || !intern->idb_ptr || !intern->idb_meth)
	return false;
    if (intern->idb_meth->magic != RT_FUNCTAB_MAGIC)
	return false;

    const uint32_t expected = intern->idb_meth->ft_internal_magic;
    if (!expected)
	return true;

    return internal_payload_magic(intern) == expected;
}

struct owned_leaf_internal {
    struct rt_db_internal local;
    struct rt_db_internal *intern;
    bool ownsLocal;

    owned_leaf_internal(void) : intern(NULL), ownsLocal(false)
    {
	RT_DB_INTERNAL_INIT(&local);
    }

    ~owned_leaf_internal(void)
    {
	if (ownsLocal)
	    rt_db_free_internal(&local);
    }
};

static struct rt_db_internal *
import_walk_leaf_internal(struct db_tree_state *tsp,
			  struct directory *dp,
			  struct owned_leaf_internal *handle)
{
    if (!handle || !tsp || !tsp->ts_dbip || !dp)
	return NULL;

    if (rt_db_get_internal(&handle->local, dp, tsp->ts_dbip, NULL) < 0)
	return NULL;

    /* db_walk_tree_leaf_instances deliberately has not imported this leaf.
     * This is the one owned primitive import for the occurrence.  Record
     * ownership before validation so even a malformed payload is released. */
    handle->ownsLocal = true;
    if (!internal_payload_magic_valid(&handle->local))
	return NULL;

    handle->intern = &handle->local;
    return handle->intern;
}

static void
set_leaf_import_diagnostic(struct realize_walk_data *data,
			   const struct db_full_path *pathp,
			   const struct rt_db_internal *intern)
{
    if (!data)
	return;

    const char *typeLabel = primitive_type_label(intern);
    char reason[256] = {0};
    if (intern && intern->idb_ptr && intern->idb_meth &&
	intern->idb_meth->ft_internal_magic) {
	snprintf(reason, sizeof(reason),
		 "invalid primitive payload for type '%s' (magic 0x%08x, expected 0x%08x)",
		 typeLabel,
		 (unsigned int)internal_payload_magic(intern),
		 (unsigned int)intern->idb_meth->ft_internal_magic);
    } else {
	snprintf(reason, sizeof(reason),
		 "primitive import failed or returned an invalid payload for type '%s'",
		 typeLabel);
    }
    set_walk_diagnostic(data, pathp, reason);
}

static void
assign_annotation_record(SoBRLVListShape *shape,
			 const struct rt_annot_internal *annot)
{
    if (!shape || !annot)
	return;

    RT_ANNOT_CK_MAGIC(annot);

    shape->sourceType = "annotation";
    shape->geometryKind = "annotation";
    shape->annotationBasePoint = SbVec3f(
				     static_cast<float>(annot->V[X]),
				     static_cast<float>(annot->V[Y]),
				     static_cast<float>(annot->V[Z]));

    const int pointCount = (annot->vert_count > static_cast<size_t>(INT_MAX)) ?
			   INT_MAX : static_cast<int>(annot->vert_count);
    if (pointCount > 0 && annot->verts) {
	std::vector<double> points(static_cast<size_t>(pointCount) * 3);
	for (int i = 0; i < pointCount; i++) {
	    const size_t offset = static_cast<size_t>(i) * 3;
	    points[offset + 0] = annot->verts[i][X];
	    points[offset + 1] = annot->verts[i][Y];
	    points[offset + 2] = 0.0;
	}
	shape->setPreciseAnnotationPoints(points.data(), pointCount);
    } else {
	shape->setPreciseAnnotationPoints(NULL, 0);
    }

    const int segmentCount = (annot->ant.count > static_cast<size_t>(INT_MAX)) ?
			     INT_MAX : static_cast<int>(annot->ant.count);
    shape->annotationSegmentTextValid.setNum(segmentCount);
    shape->annotationSegmentKind.setNum(segmentCount);
    shape->annotationSegmentStart.setNum(segmentCount);
    shape->annotationSegmentEnd.setNum(segmentCount);
    shape->annotationTextRefPoint.setNum(segmentCount);
    shape->annotationText.setNum(segmentCount);

    for (int i = 0; i < segmentCount; i++) {
	shape->annotationSegmentTextValid.set1Value(i, FALSE);
	shape->annotationSegmentKind.set1Value(i,
					       SoBRLVListShape::ANNOTATION_SEGMENT_NONE);
	shape->annotationSegmentStart.set1Value(i, 0);
	shape->annotationSegmentEnd.set1Value(i, 0);
	shape->annotationTextRefPoint.set1Value(i, 0);
	shape->annotationText.set1Value(i, "");

	const uint32_t *magic = annot->ant.segments ?
				static_cast<const uint32_t *>(annot->ant.segments[i]) : NULL;
	if (!magic)
	    continue;

	switch (*magic) {
	    case CURVE_LSEG_MAGIC: {
		const struct line_seg *lsg =
			reinterpret_cast<const struct line_seg *>(magic);
		shape->annotationSegmentKind.set1Value(i,
						       SoBRLVListShape::ANNOTATION_SEGMENT_LINE);
		shape->annotationSegmentStart.set1Value(i, lsg->start);
		shape->annotationSegmentEnd.set1Value(i, lsg->end);
		break;
	    }
	    case ANN_TSEG_MAGIC: {
		const struct txt_seg *tsg =
			reinterpret_cast<const struct txt_seg *>(magic);
		const char *label = BU_VLS_IS_INITIALIZED(&tsg->label) ?
				    bu_vls_cstr(&tsg->label) : "";
		shape->annotationSegmentKind.set1Value(i,
						       SoBRLVListShape::ANNOTATION_SEGMENT_TEXT);
		shape->annotationTextRefPoint.set1Value(i, tsg->ref_pt);
		shape->annotationText.set1Value(i, label ? label : "");
		shape->annotationSegmentTextValid.set1Value(i,
			(label && label[0]) ? TRUE : FALSE);
		break;
	    }
	    default:
		break;
	}
    }
}

static void
free_annotation_record_copy(struct rt_annot_internal *annot)
{
    if (!annot)
	return;

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    intern.idb_major_type = DB5_MAJORTYPE_BRLCAD;
    intern.idb_type = ID_ANNOT;
    intern.idb_meth = &OBJ[ID_ANNOT];
    intern.idb_ptr = annot;
    rt_db_free_internal(&intern);
}

template <typename ShapeT>
static void
apply_source_database_metadata(ShapeT *shape,
			       const SoBRLDatabaseSource *source)
{
    if (!shape || !source || !source->databaseMetadataValid.getValue())
	return;

    shape->regionId = source->databaseRegionId.getValue();
    shape->airCode = source->databaseAirCode.getValue();
    shape->materialId = source->databaseMaterialId.getValue();
    shape->los = source->databaseLos.getValue();

    shape->materialColorValid =
	source->databaseMaterialColorValid.getValue();
    if (source->databaseMaterialColorValid.getValue()) {
	shape->materialColor = source->databaseMaterialColor.getValue();
    } else if (!source->materialColorValid.getValue()) {
	shape->materialColor = SbColor(1.0f, 1.0f, 1.0f);
    }
    shape->materialShader = source->databaseMaterialShader.getValue();

    if (source->materialColorValid.getValue() &&
	(source->materialPolicy.getValue() !=
	 SoBRLDatabaseSource::MATERIAL_DATABASE ||
	 !shape->materialColorValid.getValue())) {
	shape->materialColorValid = TRUE;
	shape->materialColor = source->materialColor.getValue();
	shape->materialRevision = source->materialRevision.getValue();
    }
}

template <typename ShapeT>
static void
assign_realized_identity(ShapeT *shape,
			 const struct db_tree_state *tsp,
			 const char *path,
			 const char *sourceName,
			 const char *sourceType,
			 uint32_t sourceId,
			 const SoBRLDatabaseSource *source,
			 BObolMaterialColorSweep *materialSweep = nullptr)
{
    if (!shape)
	return;

    shape->sourcePath = path ? path : "";
    shape->sourceName = sourceName ? sourceName : "";
    shape->sourceType = sourceType ? sourceType : "";
    shape->sourceId = sourceId;
    if (source && source->displayName.getValue().getLength() > 0)
	shape->displayName = source->displayName.getValue();
    else
	shape->displayName = sourceName ? sourceName : (path ? path : "");
    shape->geometryName = sourceName ? sourceName : "";
    shape->sourceIdentity = source_record_identity(source, path);
    shape->cacheIdentity = record_identity_with_revision(
			       shape->sourceIdentity.getValue().getString(), sourceId);
    shape->databaseIntent = TRUE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = FALSE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = FALSE;
    if (shape->databaseIntent.getValue() ||
	!shape->nonDatabaseSource.getValue())
	shape->drawMode = source_record_draw_mode(source);
    shape->hiddenLine = (source_record_draw_mode(source) ==
			 BOBOL_LOD_DRAW_HIDDEN_LINE) ? TRUE : FALSE;
    shape->recordRole = "database";
    shape->geometryKind = "";
    if (source) {
	shape->visible = source->visible.getValue();
	shape->selected = source->selected.getValue();
	shape->highlighted = source->highlighted.getValue();
	shape->lineStyle = source->lineStyle.getValue();
	shape->lineWidth = source->lineWidth.getValue();
	shape->transparency = source->transparency.getValue();
	shape->colorOverride = source->colorOverride.getValue();
	shape->color = source->color.getValue();
	shape->materialColorValid = source->materialColorValid.getValue();
	shape->materialColor = source->materialColor.getValue();
	shape->materialRevision = source->materialRevision.getValue();
    }

    if (!tsp) {
	shape->regionId = 0;
	shape->airCode = 0;
	shape->materialId = 0;
	shape->los = 0;
	if (!source || !source->materialColorValid.getValue()) {
	    shape->materialColorValid = FALSE;
	    shape->materialColor = SbColor(1.0f, 1.0f, 1.0f);
	}
	shape->materialShader = "";
	apply_source_database_metadata(shape, source);
	return;
    }

    shape->regionId = tsp->ts_regionid;
    shape->airCode = tsp->ts_aircode;
    shape->materialId = tsp->ts_gmater;
    shape->los = tsp->ts_los;
    shape->materialColorValid = tsp->ts_mater.ma_color_valid ? TRUE : FALSE;
    shape->materialColor = SbColor(
			       static_cast<float>(tsp->ts_mater.ma_color[0]),
			       static_cast<float>(tsp->ts_mater.ma_color[1]),
			       static_cast<float>(tsp->ts_mater.ma_color[2]));
    if (source && source->materialColorValid.getValue() &&
	(source->materialPolicy.getValue() !=
	 SoBRLDatabaseSource::MATERIAL_DATABASE ||
	 !shape->materialColorValid.getValue())) {
	shape->materialColorValid = TRUE;
	shape->materialColor = source->materialColor.getValue();
    }
    shape->materialShader = tsp->ts_mater.ma_shader ? tsp->ts_mater.ma_shader : "";
    apply_source_database_metadata(shape, source);
    if (source && source->materialPolicy.getValue() == SoBRLDatabaseSource::MATERIAL_DATABASE) {
	BObolMaterialPathState resolved;
	SbColor pathColor;
	const bool haveColor = materialSweep ? materialSweep->resolve(path, resolved) :
	    bobol_database_source_path_material_color(source->getDatabase(), path, pathColor);
	if (haveColor) {
	    shape->materialColorValid = TRUE;
	    shape->materialColor = materialSweep ? resolved.color : pathColor;
	}
    }
}

template <typename ShapeT>
static void
assign_shared_geometry_identity(ShapeT *shape,
				const char *sourceName,
				const char *sourceType,
				uint32_t sourceId,
				const char *geometryKind)
{
    if (!shape)
	return;

    const char *name = sourceName ? sourceName : "";
    shape->sourcePath = name;
    shape->sourceName = name;
    shape->sourceType = sourceType ? sourceType : "";
    shape->sourceId = sourceId;
    shape->displayName = name;
    shape->geometryName = name;
    shape->sourceIdentity = name;
    shape->cacheIdentity = record_identity_with_revision(name, sourceId);
    shape->databaseIntent = TRUE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = FALSE;
    shape->sharedSource = TRUE;
    shape->nonDatabaseSource = FALSE;
    shape->recordRole = "shared-geometry";
    if (geometryKind && geometryKind[0])
	shape->geometryKind = geometryKind;
}

static BObolRealizedShapeSummary
compact_occurrence_summary(const SoBRLDatabaseSource *source,
	const char *path, const char *sourceName, const char *sourceType,
	const char *geometryKind, uint32_t sourceId, int shapeKind)
{
    BObolRealizedShapeSummary summary;
    summary.valid = TRUE;
    summary.shapeKind = shapeKind;
    summary.path = path ? path : "";
    summary.sourceName = sourceName ? sourceName : "";
    summary.sourceType = sourceType ? sourceType : "";
    summary.sourceId = sourceId;
    summary.displayName = source &&
	source->displayName.getValue().getLength() > 0 ?
	source->displayName.getValue() : summary.sourceName;
    summary.geometryName = summary.sourceName;
    summary.sourceIdentity = source_record_identity(source, path);
    summary.cacheIdentity = record_identity_with_revision(
	summary.sourceIdentity.getString(), sourceId);
    summary.databaseIntent = TRUE;
    summary.localSource = FALSE;
    summary.sharedSource = FALSE;
    summary.nonDatabaseSource = FALSE;
    summary.drawMode = source_record_draw_mode(source);
    summary.recordRole = "database";
    summary.geometryKind = geometryKind ? geometryKind : "";
    summary.visible = source ? source->visible.getValue() : TRUE;
    summary.selectable = TRUE;
    summary.selected = source ? source->selected.getValue() : FALSE;
    summary.highlighted = source ? source->highlighted.getValue() : FALSE;
    summary.hiddenLine = summary.drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    summary.lineStyle = source ? source->lineStyle.getValue() : 0;
    summary.lineWidth = source ? source->lineWidth.getValue() : 0;
    summary.transparency = source ? source->transparency.getValue() : 0.0f;
    summary.colorOverride = source ? source->colorOverride.getValue() : FALSE;
    summary.color = source ? source->color.getValue() :
	SbColor(1.0f, 1.0f, 1.0f);
    summary.materialColorValid = FALSE;
    summary.materialColor = SbColor(1.0f, 1.0f, 1.0f);
    summary.materialRevision = source ?
	source->materialRevision.getValue() : 0;
    if (source && source->databaseMetadataValid.getValue()) {
	summary.regionId = source->databaseRegionId.getValue();
	summary.airCode = source->databaseAirCode.getValue();
	summary.materialId = source->databaseMaterialId.getValue();
	summary.los = source->databaseLos.getValue();
	summary.materialColorValid =
	    source->databaseMaterialColorValid.getValue();
	if (summary.materialColorValid)
	    summary.materialColor = source->databaseMaterialColor.getValue();
	summary.materialShader = source->databaseMaterialShader.getValue();
    }
    if (source) {
	summary.ownerSourcePath = source->path.getValue();
	summary.ownerSourceInstanceKey = source_effective_instance_key(source);
	summary.ownerDrawMode = source->drawMode.getValue();
	summary.ownerSourceRevision = source->sourceRevision.getValue();
	summary.ownerInputsRevision = source->inputsRevision.getValue();
	summary.ownerViewRevision = source->viewRevision.getValue();
    }
    return summary;
}

static BObolRealizedShapeSummary
compact_occurrence_tree_summary(const SoBRLDatabaseSource *source,
	const struct db_tree_state *tsp, const struct db_full_path *fullPath,
	const char *path,
	const char *sourceName, const char *sourceType,
	const char *geometryKind, uint32_t sourceId, int shapeKind,
	BObolMaterialColorSweep *materialSweep)
{
    BObolRealizedShapeSummary summary = compact_occurrence_summary(source,
	path, sourceName, sourceType, geometryKind, sourceId, shapeKind);
    if (!tsp)
	return summary;
    summary.regionId = tsp->ts_regionid;
    summary.airCode = tsp->ts_aircode;
    summary.materialId = tsp->ts_gmater;
    summary.los = tsp->ts_los;
    summary.materialColorValid = tsp->ts_mater.ma_color_valid ? TRUE : FALSE;
    summary.materialColor = SbColor(
	static_cast<float>(tsp->ts_mater.ma_color[0]),
	static_cast<float>(tsp->ts_mater.ma_color[1]),
	static_cast<float>(tsp->ts_mater.ma_color[2]));
    summary.materialShader = tsp->ts_mater.ma_shader ?
	tsp->ts_mater.ma_shader : "";
    /* db_tree_state does not always contain the same effective color as
     * db_full_path_color (region-table fallback and some inherited BREP
     * colors are notable cases).  The prefix-cached sweep implements those
     * full-path rules without re-importing every combination for every leaf,
     * and retains the database color even while a display policy overrides it. */
    if (source) {
	BObolMaterialPathState resolved;
	SbColor databaseColor;
	const bool haveColor = materialSweep ?
	    (fullPath ? materialSweep->resolve(fullPath, resolved) :
		materialSweep->resolve(path, resolved)) :
	    bobol_database_source_path_material_color(source->getDatabase(),
		path, databaseColor);
	if (haveColor) {
	    summary.materialColorValid = TRUE;
	    summary.materialColor = materialSweep ? resolved.color :
		databaseColor;
	    if (materialSweep) {
		summary.regionId = resolved.regionId;
		summary.airCode = resolved.airCode;
		summary.materialId = resolved.materialId;
		summary.los = resolved.los;
		summary.materialShader = resolved.shader.c_str();
	    }
	}
    }
    /* Compact entries represent individual database occurrences.  Source
     * metadata belongs to the aggregate and must not replace the occurrence's
     * full-path material state resolved above. */
    return summary;
}

static void
compact_source_mesh_request_sync(BObolSourceMeshRequest &request,
	const BObolRealizedShapeSummary &summary)
{
    request.path = summary.path;
    request.sourceName = summary.sourceName;
    request.sourceType = summary.sourceType;
    request.sourceId = summary.sourceId;
    request.displayName = summary.displayName;
    request.geometryName = summary.geometryName;
    request.cacheIdentity = summary.cacheIdentity;
    request.sourceIdentity = summary.sourceIdentity;
    request.ownerSourceInstanceKey = summary.ownerSourceInstanceKey;
    request.databaseIntent = summary.databaseIntent;
    request.overlayIntent = summary.overlayIntent;
    request.hudIntent = summary.hudIntent;
    request.localSource = summary.localSource;
    request.sharedSource = summary.sharedSource;
    request.nonDatabaseSource = summary.nonDatabaseSource;
    request.drawMode = summary.drawMode;
    request.recordRole = summary.recordRole;
    request.geometryKind = summary.geometryKind;
    request.regionId = summary.regionId;
    request.airCode = summary.airCode;
    request.materialId = summary.materialId;
    request.los = summary.los;
    request.materialColorValid = summary.materialColorValid;
    request.materialColor = summary.materialColor;
    request.materialShader = summary.materialShader;
    request.selected = summary.selected;
    request.highlighted = summary.highlighted;
    request.ghosted = summary.ghosted;
    request.hiddenLine = summary.hiddenLine;
    request.editEmphasis = summary.editEmphasis;
    request.editIntentId = summary.editIntentId;
    request.editIntentRole = summary.editIntentRole;
    request.lodPolicy = summary.lodPolicy;
    request.colorOverride = summary.colorOverride;
    request.color = summary.color;
    request.transparency = summary.transparency;
}

static void
compact_summary_lod_from_source_mesh_request(
	BObolRealizedShapeSummary &summary,
	const BObolSourceMeshRequest &request)
{
    summary.lodPolicy = request.lodPolicy;
    summary.lodAvailable = request.lodAvailable ? TRUE : FALSE;
    summary.lodActiveCut = request.lodActiveCut;
    summary.lodFaceCount = request.lodFaceCount;
    summary.lodPointCount = request.lodPointCount;
    summary.lodOriginalPointCount = request.lodOriginalPointCount;
    summary.lodNormalCount = request.lodNormalCount;
    summary.lodHasSnappedPoints = request.lodHasSnappedPoints ? TRUE : FALSE;
    summary.lodHasNormals = request.lodHasNormals ? TRUE : FALSE;
    summary.lodBoundsMin = request.lodBoundsMin;
    summary.lodBoundsMax = request.lodBoundsMax;
}

template <typename ShapeT>
static void
sync_shape_owner_state(ShapeT *shape, const SoBRLDatabaseSource *source)
{
    if (!shape || !source)
	return;

    database_source_assign_string(shape->ownerSourcePath,
				  source->path.getValue());
    database_source_assign_string(shape->ownerSourceInstanceKey,
				  source_effective_instance_key(source));
    shape->ownerSourceRevision = source->sourceRevision.getValue();
    shape->ownerInputsRevision = source->inputsRevision.getValue();
    shape->ownerViewRevision = source->viewRevision.getValue();
    shape->ownerRealizedRevision = source->realizedRevision.getValue();
    shape->ownerRealizedSourceRevision =
	source->realizedSourceRevision.getValue();
    shape->ownerRealizedInputsRevision =
	source->realizedInputsRevision.getValue();
    shape->ownerRealizedViewRevision = source->realizedViewRevision.getValue();
    shape->ownerRealizationStatus = source->realizationStatus.getValue();
    database_source_assign_string(shape->ownerRealizationDiagnostic,
				  source->realizationDiagnostic.getValue());
    database_source_assign_string(shape->ownerRealizationIdentity,
				  source->realizationIdentity.getValue());
    shape->ownerSourceStale = source->stale.getValue();
    shape->ownerStaleReason = source->staleReason.getValue();
}

template <typename ShapeT>
static void
sync_shape_placement_state(ShapeT *shape, const SoBRLDatabaseSource *source)
{
    if (!shape || !source)
	return;

    shape->drawMatrixValid = source->drawMatrixValid.getValue();
    shape->drawMatrix = source->drawMatrix.getValue();
    shape->drawCenterValid = source->drawCenterValid.getValue();
    shape->drawCenter = source->drawCenter.getValue();
    shape->drawSizeValid = source->drawSizeValid.getValue();
    shape->drawSize = source->drawSize.getValue();
}

template <typename ShapeT>
static void
sync_shape_display_state(ShapeT *shape, const SoBRLDatabaseSource *source)
{
    if (!shape || !source)
	return;

    if (!shape->databaseIntent.getValue() &&
	shape->nonDatabaseSource.getValue())
	return;

    shape->drawMode = source_record_draw_mode(source);
    shape->hiddenLine = (shape->drawMode.getValue() ==
			 BOBOL_LOD_DRAW_HIDDEN_LINE) ? TRUE : FALSE;
    shape->visible = source->visible.getValue();
    shape->selected = source->selected.getValue();
    shape->highlighted = source->highlighted.getValue();
    shape->lineStyle = source->lineStyle.getValue();
    shape->lineWidth = source->lineWidth.getValue();
    shape->transparency = source->transparency.getValue();
    shape->colorOverride = source->colorOverride.getValue();
    shape->color = source->color.getValue();
    shape->selectedColor = source->selectedColor.getValue();
    shape->highlightedColor = source->highlightedColor.getValue();
    shape->ghostedColor = source->ghostedColor.getValue();
    if (source->materialColorValid.getValue()) {
	shape->materialColorValid = TRUE;
	shape->materialColor = source->materialColor.getValue();
	shape->materialRevision = source->materialRevision.getValue();
    }
    apply_source_database_metadata(shape, source);
}

template <typename ShapeT>
static void
sync_shape_display_name(ShapeT *shape, const SoBRLDatabaseSource *source)
{
    if (!shape || !source)
	return;

    if (!shape->databaseIntent.getValue() &&
	shape->nonDatabaseSource.getValue())
	return;

    const SbString &nextDisplayName = source->displayName.getValue().getLength() ?
	source->displayName.getValue() : (shape->sourceName.getValue().getLength() ?
	shape->sourceName.getValue() : shape->sourcePath.getValue());
    database_source_assign_string(shape->displayName,
				  nextDisplayName);
}

template <typename Shape>
static void
sync_shape_source_state(Shape *shape, const SoBRLDatabaseSource *source)
{
    sync_shape_owner_state(shape, source);
    sync_shape_placement_state(shape, source);
    sync_shape_display_state(shape, source);
    sync_shape_display_name(shape, source);
}

static void
sync_realized_shape_owner_state_in_node(SoNode *node,
					const SoBRLDatabaseSource *source)
{
    if (!node || !source || (node != source && node->isOfType(SoBRLDatabaseSource::getClassTypeId())))
	return;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	sync_shape_source_state(shape, source);
	return;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	sync_shape_source_state(shape, source);
	return;
    }
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    sync_realized_shape_owner_state_in_node(group->getChild(i), source);
    }
}


static union tree *
    make_nop_tree(void)
{
    union tree *tp = NULL;
    BU_GET(tp, union tree);
    RT_TREE_INIT(tp);
    tp->tr_op = OP_NOP;
    return tp;
}

static void
material_object_add_properties(SoBRLMaterialObject *object,
			       const char *group,
			       const struct bu_attribute_value_set *properties)
{
    const struct bu_attribute_value_pair *avpp = NULL;

    if (!object || !properties)
	return;

    for (BU_AVS_FOR(avpp, properties))
	object->addProperty(group, avpp->name, avpp->value);
}

static SoBRLMaterialObject *
material_object_from_internal(struct rt_material_internal *material)
{
    if (!material)
	return NULL;
    RT_CHECK_MATERIAL(material);

    SoBRLMaterialObject *object = new SoBRLMaterialObject;
    object->materialName = bu_vls_cstr(&material->name);
    object->parentName = bu_vls_cstr(&material->parent);
    object->materialSource = bu_vls_cstr(&material->source);
    material_object_add_properties(object, "physical",
				   &material->physicalProperties);
    material_object_add_properties(object, "mechanical",
				   &material->mechanicalProperties);
    material_object_add_properties(object, "optical",
				   &material->opticalProperties);
    material_object_add_properties(object, "thermal",
				   &material->thermalProperties);
    return object;
}

static void
assign_material_identity(SoBRLMaterialObject *object,
			 const char *path,
			 const char *sourceName,
			 const char *sourceType,
			 uint32_t sourceId)
{
    if (!object)
	return;

    object->sourcePath = path ? path : "";
    object->sourceName = sourceName ? sourceName : "";
    object->sourceType = sourceType ? sourceType : "";
    object->sourceId = sourceId;
}

static SoBRLVListShape *
vlist_from_plot_internal(struct rt_db_internal *intern,
			 const SoBRLDatabaseSource *source,
			 const struct bg_tess_tol *plotTtol = NULL,
			 const struct bn_tol *plotTol = NULL)
{
    if (!internal_payload_magic_valid(intern))
	return NULL;
    if (!intern->idb_meth || !intern->idb_meth->ft_plot)
	return NULL;

    struct rt_annot_internal *annotation = NULL;
    if (primitive_is_annotation(intern->idb_type, primitive_type_label(intern))) {
	struct rt_annot_internal *annot =
		static_cast<struct rt_annot_internal *>(intern->idb_ptr);
	RT_ANNOT_CK_MAGIC(annot);
	annotation = rt_copy_annot(annot);
    }

    std::unique_ptr<rt_annot_internal, decltype(&free_annotation_record_copy)> annotationOwner(
	annotation, free_annotation_record_copy);
    struct bg_tess_tol ttol = plotTtol ? *plotTtol : source_tess_tol(source);
    struct bn_tol tol = BN_TOL_INIT_TOL;
    if (plotTol)
	tol = *plotTol;
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    const int ret = plot_internal_to_vlist_geometry(points, commands,
	intern, &ttol, &tol);
    if (ret < 0) {
	return NULL;
    }

    if (points.empty() || points.size() != commands.size()) {
	return NULL;
    }

    SoBRLVListShape *shape = new SoBRLVListShape;
    SbModernUtils::SoNodeRef prepared(shape);
    shape->setLineSet(points.data(), commands.data(),
		      static_cast<int>(points.size()));
    assign_annotation_record(shape, annotation);
    prepared.release()->unrefNoDelete();
    return shape;
}

static int
cad_wire_part_geometry_from_plot_internal(struct rt_db_internal *intern,
	const SoBRLDatabaseSource *source, Obol::PartGeometryBuilder &geometry)
{
    if (!internal_payload_magic_valid(intern) || !intern->idb_meth ||
	!intern->idb_meth->ft_plot)
	return 0;

    struct bg_tess_tol ttol = source_tess_tol(source);
    struct bn_tol tol = BN_TOL_INIT_TOL;
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    if (plot_internal_to_vlist_geometry(points, commands, intern,
	    &ttol, &tol) < 0)
	return 0;
    return cad_wire_part_geometry_from_line_set(points, commands, geometry);
}

static double
cad_wire_point_segment_distance_squared(const SbVec3f &point,
	const SbVec3f &start, const SbVec3f &end)
{
    const SbVec3f delta = end - start;
    const double lengthSquared = static_cast<double>(delta.sqrLength());
    if (!(lengthSquared > 0.0))
	return static_cast<double>((point - start).sqrLength());
    double t = static_cast<double>((point - start).dot(delta)) /
	lengthSquared;
    t = std::max(0.0, std::min(1.0, t));
    const SbVec3f nearest = start + delta * static_cast<float>(t);
    return static_cast<double>((point - nearest).sqrLength());
}

static std::vector<size_t>
cad_wire_simplify_polyline(const std::vector<SbVec3f> &points,
	size_t first, size_t count, double tolerance)
{
    std::vector<size_t> result;
    if (count < 3 || !(tolerance > 0.0)) {
	result.resize(count);
	std::iota(result.begin(), result.end(), first);
	return result;
    }
    std::vector<unsigned char> keep(count, 0);
    keep.front() = 1;
    keep.back() = 1;
    std::vector<std::pair<size_t, size_t>> work;
    work.emplace_back(first, first + count - 1);
    const double toleranceSquared = tolerance * tolerance;
    while (!work.empty()) {
	const std::pair<size_t, size_t> range = work.back();
	work.pop_back();
	double maximum = toleranceSquared;
	size_t split = range.first;
	for (size_t i = range.first + 1; i < range.second; ++i) {
	    const double distance = cad_wire_point_segment_distance_squared(
		points[i], points[range.first], points[range.second]);
	    if (distance > maximum) {
		maximum = distance;
		split = i;
	    }
	}
	if (split != range.first) {
	    keep[split - first] = 1;
	    if (split > range.first + 1)
		work.emplace_back(range.first, split);
	    if (range.second > split + 1)
		work.emplace_back(split, range.second);
	}
    }
    result.reserve(count);
    for (size_t i = 0; i < keep.size(); ++i)
	if (keep[i])
	    result.push_back(first + i);
    return result;
}

struct cad_wire_curve_range {
    size_t first;
    size_t count;
};

static bool
cad_wire_level_equal(const std::vector<SbVec3f> &leftPoints,
	const std::vector<uint32_t> &leftIds,
	const std::vector<SbVec3f> &rightPoints,
	const std::vector<uint32_t> &rightIds)
{
    if (leftPoints.size() != rightPoints.size() || leftIds != rightIds)
	return false;
    for (size_t i = 0; i < leftPoints.size(); ++i)
	if (leftPoints[i] != rightPoints[i])
	    return false;
    return true;
}

static int
cad_progressive_wire_part_geometry_from_provider(
    struct rt_db_internal *intern, const SoBRLDatabaseSource *source,
    Obol::PartGeometryBuilder &geometry)
{
    if (!source_mesh_lod_active(source) ||
	!internal_payload_magic_valid(intern) || !intern->idb_meth ||
	!intern->idb_meth->ft_wireframe_line_set)
	return 0;

    struct rt_primitive_lod_realization realization;
    memset(&realization, 0, sizeof(realization));
    const struct bg_tess_tol ttol = source_tess_tol(source);
    const struct bn_tol tol = BN_TOL_INIT_TOL;
    const int ret = intern->idb_meth->ft_wireframe_line_set(
	&realization, intern, &ttol, &tol);
    if (ret < 0 || !realization.has_line_set ||
	!realization.line_points || !realization.line_count) {
	primitive_realization_line_set_free(&realization);
	return 0;
    }

    size_t curve_count = 0;
    for (size_t i = 0; i < realization.line_count; ++i) {
	const int command = realization.line_commands ?
	    realization.line_commands[i] : RT_PRIMITIVE_LINE_DRAW;
	if (command == RT_PRIMITIVE_LINE_MOVE)
	    curve_count++;
    }
    if (!curve_count && realization.line_count)
	curve_count = 1;

    std::vector<SbVec3f> points;
    std::vector<cad_wire_curve_range> curves;
    points.reserve(realization.line_count);
    curves.reserve(curve_count);
    for (size_t i = 0; i < realization.line_count; ++i) {
	const int command = realization.line_commands ?
	    realization.line_commands[i] : RT_PRIMITIVE_LINE_DRAW;
	const SbVec3f point(
	    static_cast<float>(realization.line_points[i][X]),
	    static_cast<float>(realization.line_points[i][Y]),
	    static_cast<float>(realization.line_points[i][Z]));
	if (!std::isfinite(point[0]) || !std::isfinite(point[1]) ||
	    !std::isfinite(point[2])) {
	    if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
		bu_log("[obol-timing] progressive wire provider returned "
		       "non-finite point %zu: %.17g %.17g %.17g\n", i,
		       static_cast<double>(realization.line_points[i][X]),
		       static_cast<double>(realization.line_points[i][Y]),
		       static_cast<double>(realization.line_points[i][Z]));
	    primitive_realization_line_set_free(&realization);
	    return 0;
	}
	if (command == RT_PRIMITIVE_LINE_MOVE || curves.empty())
	    curves.push_back({points.size(), 0});
	const cad_wire_curve_range &curve = curves.back();
	if (!curve.count || points[curve.first + curve.count - 1] != point) {
	    points.push_back(point);
	    curves.back().count++;
	}
    }
    primitive_realization_line_set_free(&realization);
    curves.erase(std::remove_if(curves.begin(), curves.end(),
	[](const cad_wire_curve_range &curve) { return curve.count < 2; }),
	curves.end());
    if (curves.empty())
	return 0;

    Obol::WireRep wire;
    wire.bounds.makeEmpty();
    for (const cad_wire_curve_range &curve : curves)
	for (size_t i = curve.first; i < curve.first + curve.count; ++i)
	    wire.bounds.extendBy(points[i]);
    if (wire.bounds.isEmpty())
	return 0;

    std::vector<SbVec3f> previousPoints;
    std::vector<uint32_t> previousIds;
    uint32_t previousFirst = 0;
    uint32_t previousCount = 0;
    wire.progressiveCuts.resize(BOBOL_MESH_LOD_CUT_COUNT_MAX);
    const int final_level = BOBOL_MESH_LOD_CUT_COUNT_MAX - 1;
    for (int level = 0; level < BOBOL_MESH_LOD_CUT_COUNT_MAX; ++level) {
	std::vector<SbVec3f> levelPoints;
	std::vector<uint32_t> levelIds;
	bool level_exact = true;
	for (size_t curveIndex = 0; curveIndex < curves.size(); ++curveIndex) {
	    const cad_wire_curve_range &curve = curves[curveIndex];
	    SbBox3f curveBounds;
	    curveBounds.makeEmpty();
	    for (size_t i = curve.first; i < curve.first + curve.count; ++i)
		curveBounds.extendBy(points[i]);
	    const SbVec3f extent = curveBounds.getMax() - curveBounds.getMin();
	    const double diagonal = sqrt(static_cast<double>(extent.sqrLength()));
	    const double tolerance = level == final_level ? 0.0 :
		ldexp(diagonal, -level - 2);
	    const std::vector<size_t> selected =
		cad_wire_simplify_polyline(points, curve.first, curve.count,
		    tolerance);
	    if (selected.size() != curve.count)
		level_exact = false;
	    else {
		for (size_t i = 0; i < selected.size(); ++i) {
		    if (selected[i] != curve.first + i) {
			level_exact = false;
			break;
		    }
		}
	    }
	    for (size_t i = 1; i < selected.size(); ++i) {
		levelPoints.push_back(points[selected[i - 1]]);
		levelPoints.push_back(points[selected[i]]);
		levelIds.push_back(static_cast<uint32_t>(
		    std::min<size_t>(curveIndex, UINT32_MAX)));
	    }
	}
	if (levelPoints.size() / 2 > UINT32_MAX ||
	    wire.segmentCount() > UINT32_MAX - levelPoints.size() / 2)
	    return 0;
	wire.progressiveCuts[level].maximumNormalizedError =
	    (level == final_level || level_exact) ?
	    0.0f : std::ldexp(1.0f, -static_cast<int>(level) - 2);
	uint32_t current_first = previousFirst;
	uint32_t current_count = previousCount;
	if (level > 0 && cad_wire_level_equal(levelPoints, levelIds,
		previousPoints, previousIds)) {
	    wire.progressiveCuts[level].segmentFirst = previousFirst;
	    wire.progressiveCuts[level].segmentCount = previousCount;
	} else {
	    current_first = static_cast<uint32_t>(wire.segmentCount());
	    current_count = static_cast<uint32_t>(levelPoints.size() / 2);
	    wire.segmentPoints.insert(wire.segmentPoints.end(),
		levelPoints.begin(), levelPoints.end());
	    wire.segmentIds.insert(wire.segmentIds.end(),
		levelIds.begin(), levelIds.end());
	    wire.progressiveCuts[level].segmentFirst = current_first;
	    wire.progressiveCuts[level].segmentCount = current_count;
	    previousPoints = std::move(levelPoints);
	    previousIds = std::move(levelIds);
	    previousFirst = current_first;
	    previousCount = current_count;
	}
	if (level_exact) {
	    for (int remaining = level + 1;
		 remaining < BOBOL_MESH_LOD_CUT_COUNT_MAX; ++remaining) {
		wire.progressiveCuts[remaining].segmentFirst = current_first;
		wire.progressiveCuts[remaining].segmentCount = current_count;
		wire.progressiveCuts[remaining].maximumNormalizedError = 0.0f;
	    }
	    break;
	}
    }
    wire.progressiveMinimumCut = 0;
    wire.progressiveResidentCut = final_level;
    wire.progressiveQuantizationMinimum = wire.bounds.getMin();
    wire.progressiveQuantizationMaximum = wire.bounds.getMax();
    geometry.wire = std::move(wire);
    return 1;
}


static SoBRLVListShape *
vlist_from_lod_realization_internal(struct rt_db_internal *intern,
				    const SoBRLDatabaseSource *source,
				    const SbBox3f &localBounds,
				    const struct bn_tol *plotTol = NULL)
{
    if (!source_csg_lod_active(source) || !internal_payload_magic_valid(intern))
	return NULL;
    if (!intern->idb_meth || !intern->idb_meth->ft_lod_realize)
	return NULL;
    if (intern->idb_type == ID_BOT || intern->idb_type == ID_MATERIAL)
	return NULL;

    struct rt_primitive_lod_realization realization;
    memset(&realization, 0, sizeof(realization));

    struct bv_view_info viewInfo;
    source_view_info(&viewInfo, source);

    struct bn_tol tol = BN_TOL_INIT_TOL;
    if (plotTol)
	tol = *plotTol;

    const fastf_t solidSize = source_lod_solid_size(source, localBounds);
    const int ret = intern->idb_meth->ft_lod_realize(&realization, intern,
	&tol, &viewInfo, solidSize);
    if (ret < 0 || !realization.has_line_set ||
	(realization.line_count && !realization.line_points)) {
	primitive_realization_line_set_free(&realization);
	return NULL;
    }

    return vlist_from_primitive_realization_line_set(&realization, "line");
}

static int
cad_wire_part_geometry_from_lod_realization_internal(
	struct rt_db_internal *intern, const SoBRLDatabaseSource *source,
	const SbBox3f &localBounds, Obol::PartGeometryBuilder &geometry,
	bool *viewDependentCsgGeometry = NULL)
{
    if (viewDependentCsgGeometry)
	*viewDependentCsgGeometry = false;
    if (cad_progressive_wire_part_geometry_from_provider(
	    intern, source, geometry))
	return 1;
    if (!source_csg_lod_active(source) || !internal_payload_magic_valid(intern)
	|| !intern->idb_meth || !intern->idb_meth->ft_lod_realize ||
	intern->idb_type == ID_BOT || intern->idb_type == ID_MATERIAL)
	return 0;

    struct rt_primitive_lod_realization realization;
    memset(&realization, 0, sizeof(realization));
    struct bv_view_info viewInfo;
    source_view_info(&viewInfo, source);
    struct bn_tol tol = BN_TOL_INIT_TOL;
    const fastf_t solidSize = source_lod_solid_size(source, localBounds);
    const int ret = intern->idb_meth->ft_lod_realize(&realization, intern,
	&tol, &viewInfo, solidSize);
    if (ret < 0 || !realization.has_line_set ||
	(realization.line_count && !realization.line_points)) {
	primitive_realization_line_set_free(&realization);
	return 0;
    }

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    points.reserve(realization.line_count);
    commands.reserve(realization.line_count);
    for (size_t i = 0; i < realization.line_count; i++) {
	const int32_t command = primitive_realization_command_to_vlist_command(
	    realization.line_commands ? realization.line_commands[i] :
	    RT_PRIMITIVE_LINE_DRAW);
	if (command < 0) {
	    primitive_realization_line_set_free(&realization);
	    return 0;
	}
	points.push_back(SbVec3f(
	    static_cast<float>(realization.line_points[i][X]),
	    static_cast<float>(realization.line_points[i][Y]),
	    static_cast<float>(realization.line_points[i][Z])));
	commands.push_back(command);
    }
    primitive_realization_line_set_free(&realization);
    const int converted = cad_wire_part_geometry_from_line_set(points,
	commands, geometry);
    if (converted && viewDependentCsgGeometry)
	*viewDependentCsgGeometry = true;
    return converted;
}

/* A wireframe BoT intentionally keeps a mesh/PoP payload so its active prefix
 * can drive both shaded and wire presentation.  Other wireframe primitives
 * require their plotted/LoD line representation; accepting a shaded cache
 * entry merely because it is resident makes output depend on which mode was
 * drawn first. */
static bool
source_cached_mesh_matches_presentation(
	const SoBRLDatabaseSource *source, const struct directory *dp)
{
    if (!source || source_record_draw_mode(source) != BOBOL_LOD_DRAW_WIRE)
	return true;
    return dp && dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT;
}

static bool
source_cached_wire_matches_mesh_presentation(
	const SoBRLDatabaseSource *source, const char *sourceType,
	const char *geometryKind)
{
    if (!source)
	return true;
    if (source_record_draw_mode(source) == BOBOL_LOD_DRAW_WIRE) {
	/*
	 * A wireframe BoT is still a view-managed triangle mesh: its active
	 * PoP prefix supplies the edges.  A persistent plotted-vlist cache may
	 * coexist with the mesh cache after a cold run, but accepting that
	 * vlist during the authoritative warm walk discards the source-mesh
	 * request.  The warm manifest initially starts PoP correctly, then final
	 * adoption replaces it with non-LoD wire geometry and strands the view
	 * at whichever prefixes happened to arrive first.
	 *
	 * Reject BoT wire caches whenever this is the mesh-role wire
	 * presentation.  The view-LoD enable bits describe whether the current
	 * controller will vary the prefix; they are not geometry identity and
	 * are deliberately unset on some detached realization workers.  Using
	 * them here made cache arbitration differ between cold and warm walks:
	 * the worker accepted a plotted BoT vlist, discarded the source request,
	 * and its authoritative handoff replaced the live PoP-backed entry.
	 *
	 * A below-threshold BoT still belongs on the mesh path: it simply keeps
	 * its complete triangle payload and derives wire edges from that payload.
	 * Other wire primitives retain their authored/evaluated line geometry.
	 */
	if (sourceType && BU_STR_EQUAL(sourceType, "bot"))
	    return false;
	return true;
    }
    if ((geometryKind && (strstr(geometryKind, "point") ||
	    BU_STR_EQUAL(geometryKind, "annotation"))) ||
	(sourceType && (BU_STR_EQUAL(sourceType, "half") ||
	    BU_STR_EQUAL(sourceType, "sketch") ||
	    BU_STR_EQUAL(sourceType, "annot") ||
	    BU_STR_EQUAL(sourceType, "annotation"))))
	return true;
    return false;
}


static union tree *
    realize_leaf(struct db_tree_state *tsp,
		 const struct db_full_path *pathp,
		 struct directory *dp,
		 void *client_data)
{
    struct realize_walk_data *data = static_cast<struct realize_walk_data *>(client_data);
    if (!data || !data->source || !pathp || !tsp || !tsp->ts_dbip ||
	!dp)
	return TREE_NULL;
    data->visited_leaves++;

    /* A streaming realization that has been cancelled stops tessellating the
     * remaining leaves promptly; the partial result is discarded by the pump. */
    if (data->stream_sink && data->stream_sink->isCancelled())
	return make_nop_tree();

    std::string walkOccurrenceIdentity;
    uint32_t duplicateOrdinal = 0;
    if (data->compact_index) {
	const std::string identity = realize_walk_instance_identity(tsp, pathp);
	if (identity.empty()) {
	    data->compact_unsupported = 1;
	    return TREE_NULL;
	}
	if (!data->compact_seen_instances.insert(identity).second)
	    return make_nop_tree();
	walkOccurrenceIdentity = realize_walk_occurrence_identity(pathp);
	duplicateOrdinal = data->compact_occurrence_counts[
	    walkOccurrenceIdentity]++;
    }

    /*
     * The compact path is the normal aggregate publication path.  Build its
     * immutable Obol geometry directly rather than allocating a Coin vlist
     * carrier and immediately converting it again below.
     */
    if (data->compact_index) {
	SbBox3f cacheBounds;
	std::string cacheKey = realize_geometry_cache_key(dp);
	source_lod_cache_key_append(cacheKey, data->source, cacheBounds,
	    dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP);
	const BObolCachedPartGeometry *cachedCad =
	    find_wire_cad_geometry_any(data->cache, cacheKey);
	std::shared_ptr<const Obol::PartGeometry> cadGeometry = cachedCad ?
	    cachedCad->geometry : std::shared_ptr<const Obol::PartGeometry>();
	bool viewDependentCsgGeometry = cachedCad ?
	    cachedCad->viewDependentCsgGeometry : false;
	bobol_performance_counter_add(cadGeometry ?
	    BOBOL_PERF_WIRE_CACHE_HITS : BOBOL_PERF_WIRE_CACHE_MISSES, 1);
	const char *typeLabel = cachedCad && !cachedCad->sourceType.empty() ?
	    cachedCad->sourceType.c_str() : (cadGeometry ? "wire" : NULL);
	const char *geometryKind = cachedCad && !cachedCad->geometryKind.empty() ?
	    cachedCad->geometryKind.c_str() : "line";
	if (!cadGeometry) {
	    owned_leaf_internal validInternal;
	    struct rt_db_internal *localIntern =
		import_walk_leaf_internal(tsp, dp, &validInternal);
	    if (!localIntern) {
		data->failed_shapes++;
		set_leaf_import_diagnostic(data, pathp,
			validInternal.ownsLocal ? &validInternal.local : NULL);
		return TREE_NULL;
	    }

	    typeLabel = primitive_type_label(localIntern);
	    SbBox3f localBounds;
	    if (source_view_lod_active(data->source))
		(void)local_bounds_from_internal(localIntern, localBounds);
	    if (localIntern->idb_type == ID_MATERIAL) {
		SoBRLMaterialObject *materialObject = material_object_from_internal(
		    static_cast<struct rt_material_internal *>(localIntern->idb_ptr));
		if (!materialObject) {
		    data->failed_shapes++;
		    set_walk_diagnostic(data, pathp,
			"material object realization failed");
		    return TREE_NULL;
		}
		char *path = db_path_to_string(pathp);
		SoSeparator *leaf = new SoSeparator;
		assign_material_identity(materialObject, path, dp->d_namep,
		    typeLabel, data->revision);
		leaf->addChild(materialObject);
		database_source_add_realized_child(data->source, leaf);
		data->realized_shapes++;
		if (path)
		    bu_free(path, "db_path_to_string");
		return make_nop_tree();
	    }

	    Obol::PartGeometryBuilder generated;
	    viewDependentCsgGeometry = false;
	    int generatedGeometry = localIntern->idb_type == ID_BOT ?
		cad_wire_part_geometry_from_bot(
		    static_cast<const struct rt_bot_internal *>(localIntern->idb_ptr),
		    generated) : cad_wire_part_geometry_from_lod_realization_internal(
			localIntern, data->source, localBounds, generated,
			&viewDependentCsgGeometry);
	    if (!generatedGeometry)
		generatedGeometry = cad_wire_part_geometry_from_plot_internal(
		    localIntern, data->source, generated);
	    if (!generatedGeometry) {
		char reason[256] = {0};
		data->failed_shapes++;
		snprintf(reason, sizeof(reason),
		    "wireframe plot produced no usable geometry for primitive type '%s'",
		    typeLabel ? typeLabel : "");
		set_walk_diagnostic(data, pathp, reason);
		return TREE_NULL;
	    }
	    if (primitive_is_annotation(localIntern->idb_type, typeLabel)) {
		typeLabel = "annotation";
		geometryKind = "annotation";
	    }
	    cadGeometry = data->cache->storeWireCadGeometry(cacheKey,
		std::move(generated), typeLabel, geometryKind,
		localBounds.isEmpty() ? NULL : &localBounds, false, NULL,
		viewDependentCsgGeometry);
	}

	char *path = db_path_to_string(pathp);
	compact_occurrence_build input;
	input.occurrence.geometry = cadGeometry;
	input.occurrence.viewDependentCsgGeometry =
	    viewDependentCsgGeometry ? TRUE : FALSE;
	input.occurrence.localTransform = mat_to_sbmatrix(tsp->ts_mat);
	input.occurrence.summary = compact_occurrence_tree_summary(
	    data->source, tsp, pathp, path, dp->d_namep,
	    geometryKind && BU_STR_EQUAL(geometryKind, "annotation") ?
	    "annotation" : typeLabel,
	    geometryKind, data->revision,
	    BObolRealizedShapeSummary::SHAPE_VLIST,
	    static_cast<BObolMaterialColorSweep *>(data->material_sweep));
	input.occurrence.occurrenceIndex =
	    data->source->occurrenceIndex.getValue();
	input.occurrence.booleanOperation =
	    data->source->booleanOperation.getValue();
	input.semantic = compact_semantic_from_summary(input.occurrence.summary);
	input.dashed = (tsp->ts_sofar & TS_SOFAR_MINUS) ? TRUE : FALSE;
	const size_t entryCount = data->compact_index->entries.size();
	compact_add_occurrence(data->source, *data->compact_index, input,
	    data->compact_ordinal, data->compact_unsupported);
	compact_apply_walk_identity(data->source, *data->compact_index,
	    entryCount, tsp, pathp, walkOccurrenceIdentity, duplicateOrdinal);
	if (data->compact_index->entries.size() > entryCount)
	    realize_walk_stream_push_current(data, input.occurrence,
		*data->compact_index, entryCount);
	SbBox3f bounds = database_source_transform_bounds(
	    compact_part_geometry_bounds(cadGeometry),
	    input.occurrence.localTransform);
	realize_walk_extend_bounds(data, bounds);
	data->realized_shapes++;
	if (path)
	    bu_free(path, "db_path_to_string");
	return make_nop_tree();
    }

    SoBRLVListShape *sharedShape = NULL;
    SbBox3f cacheBounds;
    std::string cacheKey = realize_geometry_cache_key(dp);
    source_lod_cache_key_append(cacheKey, data->source, cacheBounds,
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP);
    BObolRealizationCacheMap<SoBRLVListShape *>::iterator found =
	data->cache->sharedWireGeometry.find(cacheKey);
    if (found != data->cache->sharedWireGeometry.end())
	sharedShape = found->second;
    std::shared_ptr<const Obol::PartGeometry> cachedCadGeometry;
    const BObolCachedPartGeometry *cachedCad = NULL;
    if (data->compact_index) {
	cachedCad = find_wire_cad_geometry_any(data->cache, cacheKey);
	if (cachedCad)
	    cachedCadGeometry = cachedCad->geometry;
    }
    bobol_performance_counter_add(
	(sharedShape || cachedCadGeometry) ? BOBOL_PERF_WIRE_CACHE_HITS :
	BOBOL_PERF_WIRE_CACHE_MISSES, 1);

    const char *typeLabel = sharedShape ?
			    sharedShape->sourceType.getValue().getString() :
			    (cachedCad && !cachedCad->sourceType.empty() ?
			     cachedCad->sourceType.c_str() :
			     (cachedCadGeometry ? "wire" : NULL));
    if (!sharedShape && !cachedCadGeometry) {
	owned_leaf_internal validInternal;
	struct rt_db_internal *localIntern =
	    import_walk_leaf_internal(tsp, dp, &validInternal);
	if (!localIntern) {
	    data->failed_shapes++;
	    set_leaf_import_diagnostic(data, pathp,
		validInternal.ownsLocal ? &validInternal.local : NULL);
	    return TREE_NULL;
	}

	typeLabel = primitive_type_label(localIntern);
	SbBox3f localBounds;
	if (source_view_lod_active(data->source))
	    (void)local_bounds_from_internal(localIntern, localBounds);
	if (localIntern->idb_type == ID_MATERIAL) {
	    SoBRLMaterialObject *materialObject =
		material_object_from_internal(static_cast<struct rt_material_internal *>(localIntern->idb_ptr));
	    if (!materialObject) {
		data->failed_shapes++;
		set_walk_diagnostic(data, pathp,
				    "material object realization failed");
		return TREE_NULL;
	    }
	    char *path = db_path_to_string(pathp);
	    SoSeparator *leaf = new SoSeparator;
	    assign_material_identity(materialObject, path, dp->d_namep,
				     typeLabel, data->revision);
	    leaf->addChild(materialObject);
	    database_source_add_realized_child(data->source, leaf);
	    data->realized_shapes++;
	    if (path)
		bu_free(path, "db_path_to_string");
	    return make_nop_tree();
	}

	if (localIntern->idb_type == ID_BOT) {
	    sharedShape = vlist_from_bot_wireframe(
		static_cast<const struct rt_bot_internal *>(
		    localIntern->idb_ptr));
	} else {
	    sharedShape = vlist_from_lod_realization_internal(localIntern,
		data->source, localBounds);
	    if (!sharedShape)
		sharedShape = vlist_from_plot_internal(localIntern, data->source);
	}
	if (!sharedShape) {
	    char reason[256] = {0};
	    data->failed_shapes++;
	    snprintf(reason, sizeof(reason),
		     "wireframe plot produced no usable geometry for primitive type '%s'",
		     typeLabel);
	    set_walk_diagnostic(data, pathp, reason);
	    return TREE_NULL;
	}
	assign_shared_geometry_identity(sharedShape, dp->d_namep, typeLabel,
					data->revision, "line");
	if (primitive_is_annotation(localIntern->idb_type, typeLabel)) {
	    sharedShape->sourceType = "annotation";
	    sharedShape->geometryKind = "annotation";
	}
	if (!data->compact_index)
	    data->cache->storeWireGeometry(cacheKey, sharedShape);
	typeLabel = sharedShape->sourceType.getValue().getString();
    }

    char *path = db_path_to_string(pathp);
    if (data->compact_index) {
	const char *geometryKind = sharedShape ?
	    sharedShape->geometryKind.getValue().getString() :
	    (cachedCad && !cachedCad->geometryKind.empty() ?
	     cachedCad->geometryKind.c_str() : "line");
	std::shared_ptr<const Obol::PartGeometry> cadGeometry =
	    cachedCadGeometry;
	if (!cadGeometry && sharedShape) {
	    Obol::PartGeometryBuilder generated;
	    if (cad_vlist_part_geometry(sharedShape, generated))
		cadGeometry = data->cache->storeWireCadGeometry(cacheKey,
		    std::move(generated), typeLabel, geometryKind);
	}
	if (!cadGeometry) {
	    data->compact_unsupported = 1;
	} else {
	    compact_occurrence_build input;
	    input.occurrence.geometry = cadGeometry;
	    input.occurrence.localTransform = mat_to_sbmatrix(tsp->ts_mat);
	    input.occurrence.summary = compact_occurrence_tree_summary(
		data->source, tsp, pathp, path, dp->d_namep,
		geometryKind && BU_STR_EQUAL(geometryKind, "annotation") ?
		"annotation" : typeLabel,
		geometryKind, data->revision,
		BObolRealizedShapeSummary::SHAPE_VLIST,
		static_cast<BObolMaterialColorSweep *>(data->material_sweep));
	    input.occurrence.occurrenceIndex =
		data->source->occurrenceIndex.getValue();
	    input.occurrence.booleanOperation =
		data->source->booleanOperation.getValue();
	    input.semantic = compact_semantic_from_summary(
		input.occurrence.summary);
	    input.dashed = (tsp->ts_sofar & TS_SOFAR_MINUS) ? TRUE : FALSE;
	    const size_t entryCount = data->compact_index->entries.size();
	    compact_add_occurrence(data->source, *data->compact_index, input,
		data->compact_ordinal, data->compact_unsupported);
	    compact_apply_walk_identity(data->source, *data->compact_index,
		entryCount, tsp, pathp, walkOccurrenceIdentity,
		duplicateOrdinal);
	    if (data->compact_index->entries.size() > entryCount)
		realize_walk_stream_push_current(data, input.occurrence,
		    *data->compact_index, entryCount);
	    SbBox3f bounds = database_source_transform_bounds(
		compact_part_geometry_bounds(cadGeometry),
		input.occurrence.localTransform);
	    realize_walk_extend_bounds(data, bounds);
	}
    } else {
	BObolPerformanceTimer timer(BOBOL_PERF_REALIZED_INSTANCE_NODE_US);
	SoSeparator *leaf = realize_instance_leaf_separator(tsp);
	SoBRLVListShape *shape = new SoBRLVListShape;
	assign_realized_identity(shape, tsp, path, dp->d_namep, typeLabel,
	    data->revision, data->source,
	    static_cast<BObolMaterialColorSweep *>(data->material_sweep));
	shape->setSharedGeometry(sharedShape);
	const char *geometryKind = sharedShape->geometryKind.getValue().getString();
	shape->geometryKind = geometryKind && geometryKind[0] ? geometryKind : "line";
	if (geometryKind && BU_STR_EQUAL(geometryKind, "annotation"))
	    shape->sourceType = "annotation";
	leaf->addChild(shape);
	database_source_add_realized_child(data->source, leaf);
	if (timer.active())
	    bobol_performance_counter_add(
		BOBOL_PERF_REALIZED_INSTANCE_NODES, 1);
    }
    data->realized_shapes++;
    if (path)
	bu_free(path, "db_path_to_string");

    return make_nop_tree();
}

std::string
database_source_leaf_component(const SbString &path)
{
    const char *name = database_source_skip_leading_slash(path.getString());
    if (!name || !name[0])
	return std::string();

    const char *slash = strrchr(name, '/');
    std::string leaf((slash && slash[1]) ? slash + 1 : name);
    const size_t instanceSpecifier = leaf.find('@');
    if (instanceSpecifier != std::string::npos)
	leaf.erase(instanceSpecifier);
    return leaf;
}

bool
database_source_path_matches_frontier(const SbString &candidatePath, const char *queryPath)
{
    const char *query = database_source_skip_leading_slash(queryPath ? queryPath : "");
    if (!query[0])
	return true;
    if (!strchr(query, '/') && !strchr(query, '@') &&
	database_source_leaf_component(candidatePath) == query)
	return true;
    const char *candidate = database_source_skip_leading_slash(candidatePath.getString());
    const size_t length = strlen(query);
    if (bu_strncmp(candidate, query, length))
	return false;
    const char suffix = candidate[length];
    return suffix == '\0' || suffix == '/' || suffix == '@';
}

static std::string
database_source_full_path_string(const SbString &path)
{
    const char *pathString = path.getString();
    if (!pathString || !pathString[0])
	return std::string();
    if (pathString[0] == '/')
	return std::string(pathString);
    return std::string("/") + pathString;
}

static int
bot_lod_proxy_bounds(const struct rt_db_internal *intern,
		     uint32_t faceThreshold,
		     SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!intern || faceThreshold == 0 || intern->idb_type != ID_BOT ||
	intern->idb_ptr == NULL)
	return 0;

    const struct rt_bot_internal *bot =
	static_cast<const struct rt_bot_internal *>(intern->idb_ptr);
    if (!bot || !bot->vertices || bot->num_vertices == 0 ||
	bot->num_faces < static_cast<size_t>(faceThreshold))
	return 0;
    RT_BOT_CK_MAGIC(bot);

    for (size_t i = 0; i < bot->num_vertices; i++) {
	bounds.extendBy(SbVec3f(
	    static_cast<float>(bot->vertices[i * 3]),
	    static_cast<float>(bot->vertices[i * 3 + 1]),
	    static_cast<float>(bot->vertices[i * 3 + 2])));
    }

    return bounds.isEmpty() ? 0 : 1;
}

static SoBRLVListShape *
vlist_from_aabb_proxy_bounds(const SbBox3f &bounds)
{
    if (bounds.isEmpty())
	return NULL;

    const SbVec3f bmin = bounds.getMin();
    const SbVec3f bmax = bounds.getMax();
    const SbVec3f corners[8] = {
	SbVec3f(bmin[0], bmin[1], bmin[2]),
	SbVec3f(bmax[0], bmin[1], bmin[2]),
	SbVec3f(bmax[0], bmax[1], bmin[2]),
	SbVec3f(bmin[0], bmax[1], bmin[2]),
	SbVec3f(bmin[0], bmin[1], bmax[2]),
	SbVec3f(bmax[0], bmin[1], bmax[2]),
	SbVec3f(bmax[0], bmax[1], bmax[2]),
	SbVec3f(bmin[0], bmax[1], bmax[2])
    };
    static const int edges[12][2] = {
	{0, 1}, {1, 2}, {2, 3}, {3, 0},
	{4, 5}, {5, 6}, {6, 7}, {7, 4},
	{0, 4}, {1, 5}, {2, 6}, {3, 7}
    };

    SbVec3f points[24];
    int32_t commands[24];
    for (size_t i = 0; i < 12; i++) {
	points[i * 2] = corners[edges[i][0]];
	points[i * 2 + 1] = corners[edges[i][1]];
	commands[i * 2] = SoBRLVListShape::MOVE;
	commands[i * 2 + 1] = SoBRLVListShape::DRAW;
    }

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->setLineSet(points, commands, 24);
    shape->geometryKind = "proxy";
    return shape;
}

static SbBool
point_bbox_valid(const point_t bmin, const point_t bmax)
{
    for (int i = 0; i < 3; i++) {
	if (!isfinite(bmin[i]) || !isfinite(bmax[i]) || bmin[i] > bmax[i])
	    return FALSE;
    }
    return TRUE;
}

static SbBool
local_bounds_from_internal(struct rt_db_internal *intern, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!intern || !intern->idb_meth || !intern->idb_meth->ft_bbox)
	return FALSE;

    point_t bmin;
    point_t bmax;
    VSETALL(bmin, INFINITY);
    VSETALL(bmax, -INFINITY);
    const struct bn_tol tol = BN_TOL_INIT_TOL;
    if (intern->idb_meth->ft_bbox(intern, &bmin, &bmax, &tol) != 0 ||
	!point_bbox_valid(bmin, bmax))
	return FALSE;

    bounds = database_source_box_from_minmax(
		 SbVec3f(static_cast<float>(bmin[X]),
			 static_cast<float>(bmin[Y]),
			 static_cast<float>(bmin[Z])),
		 SbVec3f(static_cast<float>(bmax[X]),
			 static_cast<float>(bmax[Y]),
			 static_cast<float>(bmax[Z])));
    return bounds.isEmpty() ? FALSE : TRUE;
}

/* Obtain a representation-independent, Boolean-aware bound for one database
 * path without constructing raytracing regions or acceleration structures.
 * This remains detached-worker work: primitive fallback plotting and nested
 * combination imports do not belong on the GUI publication path. */
static SbBool
source_bounds_from_database_path(SoBRLDatabaseSource *source,
	const char *path, SbBox3f &bounds)
{
    bounds.makeEmpty();
    struct db_i *dbip = source ? source->getDatabase() : NULL;
    if (!dbip || !path || !path[0])
	return FALSE;

    point_t bmin;
    point_t bmax;
    struct bu_vls messages = BU_VLS_INIT_ZERO;
    const char *paths[1] = {path};
    const int ret = rt_display_bounds(&messages, dbip, 1, paths, bmin, bmax);
    bu_vls_free(&messages);
    if (ret != BRLCAD_OK || !point_bbox_valid(bmin, bmax))
	return FALSE;

    bounds = database_source_box_from_minmax(
	SbVec3f(static_cast<float>(bmin[X]), static_cast<float>(bmin[Y]),
	    static_cast<float>(bmin[Z])),
	SbVec3f(static_cast<float>(bmax[X]), static_cast<float>(bmax[Y]),
	    static_cast<float>(bmax[Z])));
    return bounds.isEmpty() ? FALSE : TRUE;
}

static void
set_source_bounds_from_local_box(SoBRLDatabaseSource *source,
				 const SbBox3f &bounds,
				 SbBool exact = TRUE)
{
    if (!source)
	return;

    if (bounds.isEmpty()) {
	source->clearSourceBounds();
	return;
    }

    (void)source->setSourceBoundsState(TRUE, bounds.getMin(),
				       bounds.getMax(), exact);
}

/* The retained payload union is an exact bound of the drawable presentation.
 * It is deliberately conservative with respect to evaluated CSG, but unlike
 * the database bound it remains available for partially drawable trees. */
static void
set_source_bounds_from_presentation(SoBRLDatabaseSource *source,
				    const SbBox3f &bounds)
{
    set_source_bounds_from_local_box(source, bounds, TRUE);
}

static int
source_has_auxiliary_children(const SoBRLDatabaseSource *source);
static int
cad_vlist_part_geometry(const SoBRLVListShape *shape,
			Obol::PartGeometryBuilder &geometry);

/* Close the representation contract for a compact registry produced by a
 * serial walk.  Parallel mesh coverage publishes these fields while work is
 * running; direct leaves and ordinary wire walks learn their exact population
 * only after it is installed. */
bool
close_compact_stream_contract(SoBRLDatabaseSource *source,
	BObolCompactOccurrenceStream *stream)
{
    if (!source || !stream || stream->isCancelled())
	return false;

    const int compactCount = source->getCompactInstanceCount();
    if (compactCount <= 0)
	return false;

    const size_t occurrenceCount = static_cast<size_t>(compactCount);
    BObolCompactSourceProfile profile;
    profile.occurrenceCount = occurrenceCount;
    std::unordered_set<std::string> assets;
    for (int occurrenceIndex = 0; occurrenceIndex < compactCount;
	 occurrenceIndex++) {
	BObolCompactOccurrence occurrence;
	if (!source->getCompactOccurrence(occurrenceIndex, occurrence))
	    return false;
	stream->recordManifestOccurrence(occurrence);
	const char *requestedAsset = occurrence.sourceMeshRequestValid ?
	    occurrence.sourceMeshRequest.meshAssetName.getString() : NULL;
	/* An empty mesh-asset identity means the occurrence's source object is
	 * canonical; this is the same fallback used by manifest persistence. */
	const char *assetName = requestedAsset && requestedAsset[0] ?
	    requestedAsset : occurrence.summary.sourceName.getString();
	if (!assetName || !assetName[0] ||
	    !assets.insert(assetName).second)
	    continue;
	struct directory *asset = db_lookup(source->getDatabase(), assetName,
	    LOOKUP_QUIET);
	const uint64_t encodedBytes = database_source_encoded_source_bytes(
	    source->getDatabase(), asset);
	profile.largestAssetBytes = std::max(profile.largestAssetBytes,
	    encodedBytes);
	profile.encodedSourceBytes = encodedBytes >
	    UINT64_MAX - profile.encodedSourceBytes ? UINT64_MAX :
	    profile.encodedSourceBytes + encodedBytes;
    }
    profile.uniqueAssetCount = assets.size();
    profile.reusedOccurrenceCount = profile.occurrenceCount >
	profile.uniqueAssetCount ? profile.occurrenceCount -
	profile.uniqueAssetCount : 0;
    profile.valid = profile.uniqueAssetCount > 0 &&
	profile.uniqueAssetCount <= profile.occurrenceCount &&
	profile.encodedSourceBytes > 0 && profile.largestAssetBytes > 0 ?
	TRUE : FALSE;

    const BObolCompactSourceProfile *certifiedProfile =
	profile.isValid(occurrenceCount) ? &profile : nullptr;
    stream->setExpectedCount(occurrenceCount);
    stream->setPreparationWorkCount(occurrenceCount);
    stream->completePreparationWork();
    if (certifiedProfile)
	stream->setSourceProfile(profile);
    (void)source->certifyCompactStream(occurrenceCount, certifiedProfile);
    return stream->sealManifest(occurrenceCount);
}

/* A direct primitive bypasses the general coverage producer, but it must
 * close the same stream contract.  In particular, the consumer needs an
 * exact population denominator and the persisted leaf manifest needs a
 * source profile.  Leaving those fields at zero made a valid one-leaf draw
 * indistinguishable from an empty/incomplete discovery result and could make
 * its first LoD request terminal while an identical region-wrapped leaf
 * worked normally. */
static void
publish_direct_compact_occurrence(SoBRLDatabaseSource *source,
	BObolCompactOccurrenceStream *stream, const struct directory *dp,
	const BObolCompactOccurrence &occurrence)
{
    if (!source || !stream || !dp || !occurrence.geometry ||
	stream->isCancelled())
	return;

    /* setCompactOccurrence canonicalizes bounds, style and source identity.
     * Publish that installed value rather than the construction input so the
     * stream journal and the live compact registry cannot disagree. */
    BObolCompactOccurrence published = occurrence;
    (void)source->getCompactOccurrence(0, published);

    (void)close_compact_stream_contract(source, stream);
    stream->push(std::move(published));
}

static int
realize_direct_leaf_wireframe(SoBRLDatabaseSource *source,
			      BObolDatabaseSourceRealizationCache *cache,
			      uint32_t revision)
{
    struct db_i *dbip = source ? source->getDatabase() : NULL;
    if (!source || !cache || !dbip)
	return 0;

    const std::string leafName =
	database_source_leaf_component(source->path.getValue());
    if (leafName.empty())
	return 0;

    struct directory *dp =
	db_lookup(dbip, leafName.c_str(), LOOKUP_QUIET);
    if (!dp || (dp->d_flags & RT_DIR_COMB))
	return 0;

    const std::string fullPath =
	database_source_full_path_string(source->path.getValue());
    SbMatrix pathMatrix;
    if (!database_source_path_matrix(dbip, source->path.getValue(),
	pathMatrix))
	return 0;
    SoBRLVListShape *sharedShape = NULL;
    SbBox3f localBounds;
    SbBool localBoundsValid = FALSE;
    std::string cacheKey = realize_geometry_cache_key(dp);
    const uint32_t wireLodThreshold = source_mesh_lod_active(source) ?
				      source->lodBotThreshold.getValue() : 0;
    if (wireLodThreshold > 0 &&
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) {
	char suffix[64] = {0};
	snprintf(suffix, sizeof(suffix), ":wire-lod-proxy:%u",
		 wireLodThreshold);
	cacheKey += suffix;
    }
    source_lod_cache_key_append(cacheKey, source, localBounds,
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP);
    BObolRealizationCacheMap<SoBRLVListShape *>::iterator found =
	cache->sharedWireGeometry.find(cacheKey);
    if (found != cache->sharedWireGeometry.end())
	sharedShape = found->second;
    bobol_performance_counter_add(sharedShape ? BOBOL_PERF_WIRE_CACHE_HITS :
	BOBOL_PERF_WIRE_CACHE_MISSES, 1);

    BObolRealizationCacheMap<SbBox3f>::const_iterator boundsFound =
	cache->sharedWireBounds.find(cacheKey);
    if (boundsFound != cache->sharedWireBounds.end() &&
	!boundsFound->second.isEmpty()) {
	localBounds = boundsFound->second;
	localBoundsValid = TRUE;
    }

    const char *typeLabel = sharedShape ?
			    sharedShape->sourceType.getValue().getString() : NULL;
    int usedLodProxy = 0;
    if (!sharedShape) {
	owned_leaf_internal validInternal;
	if (rt_db_get_internal(&validInternal.local, dp, dbip,
		NULL) < 0 ||
	    !internal_payload_magic_valid(&validInternal.local)) {
	    SbString msg;
	    msg.sprintf("%s: direct leaf wireframe internal fetch failed",
			fullPath.c_str());
	    source->realizationDiagnostic = msg;
	    if (validInternal.local.idb_ptr)
		rt_db_free_internal(&validInternal.local);
	    return -1;
	}
	validInternal.ownsLocal = true;

	typeLabel = primitive_type_label(&validInternal.local);
	if (validInternal.local.idb_type == ID_BOT ||
	    source_view_lod_active(source))
	    localBoundsValid = local_bounds_from_internal(&validInternal.local,
		localBounds);
	SbBox3f lodProxyBounds;
	usedLodProxy = bot_lod_proxy_bounds(&validInternal.local,
	    wireLodThreshold, lodProxyBounds);
	if (usedLodProxy) {
	    sharedShape = vlist_from_aabb_proxy_bounds(lodProxyBounds);
	    if (sharedShape) {
		assign_shared_geometry_identity(sharedShape, dp->d_namep,
						typeLabel, revision, "proxy");
		localBounds = lodProxyBounds;
		localBoundsValid = TRUE;
	    }
	} else if (validInternal.local.idb_type == ID_MATERIAL) {
	    SoBRLMaterialObject *materialObject =
		material_object_from_internal(
		    static_cast<struct rt_material_internal *>(
			validInternal.local.idb_ptr));
	    if (!materialObject) {
		SbString msg;
		msg.sprintf("%s: material object realization failed",
			    fullPath.c_str());
		source->realizationDiagnostic = msg;
		return -1;
	    }
	    SoSeparator *leaf = realize_matrix_leaf_separator(pathMatrix);
	    assign_material_identity(materialObject,
				     fullPath.c_str(),
				     dp->d_namep, typeLabel, revision);
	    leaf->addChild(materialObject);
	    database_source_add_realized_child(source, leaf);
	    return 1;
	}

	if (!sharedShape) {
	    if (validInternal.local.idb_type == ID_BOT) {
		sharedShape = vlist_from_bot_wireframe(
		    static_cast<const struct rt_bot_internal *>(
			validInternal.local.idb_ptr));
	    } else {
		sharedShape = vlist_from_lod_realization_internal(
		    &validInternal.local, source, localBounds);
		if (!sharedShape)
		    sharedShape = vlist_from_plot_internal(&validInternal.local,
			source);
	    }
	}
	if (!sharedShape) {
	    SbString msg;
	    msg.sprintf(
		"%s: direct leaf wireframe plot produced no usable geometry for primitive type '%s'",
		fullPath.c_str(), typeLabel ? typeLabel : "");
	    source->realizationDiagnostic = msg;
	    return -1;
	}
	if (!usedLodProxy) {
	    assign_shared_geometry_identity(sharedShape, dp->d_namep,
					    typeLabel, revision, "line");
	    if (primitive_is_annotation(validInternal.local.idb_type,
		    typeLabel)) {
		sharedShape->sourceType = "annotation";
		sharedShape->geometryKind = "annotation";
	    }
	}
	if (localBoundsValid)
	    cache->storeWireBounds(cacheKey, localBounds);
	typeLabel = sharedShape->sourceType.getValue().getString();
    }

    {
	BObolPerformanceTimer timer(BOBOL_PERF_REALIZED_INSTANCE_NODE_US);
	SoSeparator *leaf = realize_matrix_leaf_separator(pathMatrix);
	SoBRLVListShape *shape = new SoBRLVListShape;
	assign_realized_identity(shape, NULL, fullPath.c_str(),
				 dp->d_namep, typeLabel, revision, source);
	shape->setSharedGeometry(sharedShape);
	const char *geometryKind =
	    sharedShape->geometryKind.getValue().getString();
	shape->geometryKind = (geometryKind && geometryKind[0]) ?
			      geometryKind : "line";
	if (geometryKind && BU_STR_EQUAL(geometryKind, "annotation"))
	    shape->sourceType = "annotation";
	leaf->addChild(shape);
	database_source_add_realized_child(source, leaf);
	if (timer.active())
	    bobol_performance_counter_add(
		BOBOL_PERF_REALIZED_INSTANCE_NODES, 1);
    }
    if (localBoundsValid)
	set_source_bounds_from_local_box(source, localBounds);
    else
	source->clearSourceBounds();
    return 1;
}

static int
realize_direct_leaf_wireframe_compact(
    SoBRLDatabaseSource *source,
    BObolDatabaseSourceRealizationCache *cache,
    uint32_t revision,
    BObolCompactOccurrenceStream *stream)
{
    struct db_i *dbip = source ? source->getDatabase() : NULL;
    if (!source || !cache || !dbip || source_has_auxiliary_children(source))
	return 0;

    const std::string leafName =
	database_source_leaf_component(source->path.getValue());
    if (leafName.empty())
	return 0;

    struct directory *dp =
	db_lookup(dbip, leafName.c_str(), LOOKUP_QUIET);
    if (!dp || (dp->d_flags & RT_DIR_COMB))
	return 0;

    const std::string fullPath =
	database_source_full_path_string(source->path.getValue());
    SbMatrix pathMatrix;
    if (!database_source_path_matrix(dbip, source->path.getValue(),
	pathMatrix))
	return 0;
    SbBox3f localBounds;
    SbBool localBoundsValid = FALSE;
    std::string cacheKey = realize_geometry_cache_key(dp);
    const uint32_t wireLodThreshold = source_mesh_lod_active(source) ?
				      source->lodBotThreshold.getValue() : 0;
    if (wireLodThreshold > 0 &&
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) {
	char suffix[64] = {0};
	snprintf(suffix, sizeof(suffix), ":wire-lod-proxy:%u",
		 wireLodThreshold);
	cacheKey += suffix;
    }
    source_lod_cache_key_append(cacheKey, source, localBounds,
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP);
    std::shared_ptr<const Obol::PartGeometry> cadGeometry;
    const BObolCachedPartGeometry *cachedCad =
	find_wire_cad_geometry_any(cache, cacheKey);
    bool viewDependentCsgGeometry = cachedCad ?
	cachedCad->viewDependentCsgGeometry : false;
    bool lodBacked = cachedCad && cachedCad->lodBacked &&
	cachedCad->sourceMeshRequestValid;
    BObolSourceMeshRequest sourceMeshRequest;
    if (lodBacked)
	sourceMeshRequest = cachedCad->sourceMeshRequest;
    if (cachedCad)
	cadGeometry = cachedCad->geometry;
    bobol_performance_counter_add(
	(cadGeometry ? BOBOL_PERF_WIRE_CACHE_HITS :
	BOBOL_PERF_WIRE_CACHE_MISSES), 1);

    const char *typeLabel = cachedCad && !cachedCad->sourceType.empty() ?
	cachedCad->sourceType.c_str() : (cadGeometry ? "wire" : NULL);
    const char *geometryKind = cachedCad && !cachedCad->geometryKind.empty() ?
	cachedCad->geometryKind.c_str() : "line";
    if (cachedCad && !cachedCad->bounds.isEmpty()) {
	localBounds = cachedCad->bounds;
	localBoundsValid = TRUE;
    }
    if (!cadGeometry) {
	owned_leaf_internal validInternal;
	if (rt_db_get_internal(&validInternal.local, dp, dbip, NULL) < 0 ||
	    !internal_payload_magic_valid(&validInternal.local)) {
	    SbString msg;
	    msg.sprintf("%s: direct compact leaf wireframe internal fetch failed",
			fullPath.c_str());
	    source->realizationDiagnostic = msg;
	    if (validInternal.local.idb_ptr)
		rt_db_free_internal(&validInternal.local);
	    return -1;
	}
	validInternal.ownsLocal = true;

	typeLabel = primitive_type_label(&validInternal.local);
	if (validInternal.local.idb_type == ID_BOT ||
	    source_view_lod_active(source))
	    localBoundsValid = local_bounds_from_internal(&validInternal.local,
		localBounds);
	SbBox3f lodProxyBounds;
	const int usedLodProxy = bot_lod_proxy_bounds(&validInternal.local,
	    wireLodThreshold, lodProxyBounds);
	if (usedLodProxy) {
	    const struct rt_bot_internal *bot =
		static_cast<const struct rt_bot_internal *>(
		    validInternal.local.idb_ptr);
	    lodBacked = cad_source_mesh_request_from_bot(
		sourceMeshRequest, bot) != 0;
	    if (lodBacked) {
		sourceMeshRequest.meshAssetPath = fullPath.c_str();
		sourceMeshRequest.meshAssetName = dp->d_namep;
	    }
	    Obol::PartGeometryBuilder generated;
	    if (lodBacked &&
		cad_wire_part_geometry_from_aabb(lodProxyBounds, generated)) {
		geometryKind = "aabb";
		localBounds = lodProxyBounds;
		localBoundsValid = TRUE;
		cadGeometry = cache->storeWireCadGeometry(cacheKey,
		    std::move(generated), typeLabel, geometryKind, &localBounds,
		    true, &sourceMeshRequest);
	    }
	} else if (validInternal.local.idb_type == ID_MATERIAL) {
	    SoBRLMaterialObject *materialObject = material_object_from_internal(
		static_cast<struct rt_material_internal *>(
		    validInternal.local.idb_ptr));
	    if (!materialObject) {
		SbString msg;
		msg.sprintf("%s: material object realization failed",
		    fullPath.c_str());
		source->realizationDiagnostic = msg;
		return -1;
	    }
	    SoSeparator *leaf = realize_matrix_leaf_separator(pathMatrix);
	    assign_material_identity(materialObject, fullPath.c_str(),
		dp->d_namep, typeLabel, revision);
	    leaf->addChild(materialObject);
	    database_source_add_realized_child(source, leaf);
	    source->clearSourceBounds();
	    return 1;
	}

	if (!cadGeometry) {
	    Obol::PartGeometryBuilder generated;
	    viewDependentCsgGeometry = false;
	    int generatedGeometry = 0;
	    if (validInternal.local.idb_type == ID_BOT)
		generatedGeometry = cad_wire_part_geometry_from_bot(
		    static_cast<const struct rt_bot_internal *>(
			validInternal.local.idb_ptr), generated);
	    else
		generatedGeometry =
		    cad_wire_part_geometry_from_lod_realization_internal(
			&validInternal.local, source, localBounds, generated,
			&viewDependentCsgGeometry);
	    if (!generatedGeometry)
		generatedGeometry = cad_wire_part_geometry_from_plot_internal(
		    &validInternal.local, source, generated);
	    if (generatedGeometry) {
		if (primitive_is_annotation(validInternal.local.idb_type,
			typeLabel)) {
		    typeLabel = "annotation";
		    geometryKind = "annotation";
		}
		cadGeometry = cache->storeWireCadGeometry(cacheKey,
		    std::move(generated), typeLabel, geometryKind,
		    localBoundsValid ? &localBounds : NULL, false, NULL,
		    viewDependentCsgGeometry);
	    }
	}
	if (!cadGeometry) {
	    SbString msg;
	    msg.sprintf(
		"%s: direct compact leaf wireframe plot produced no usable geometry for primitive type '%s'",
		fullPath.c_str(), typeLabel ? typeLabel : "");
	    source->realizationDiagnostic = msg;
	    return -1;
	}
    }

    BObolCompactOccurrence occurrence;
    occurrence.geometry = cadGeometry;
    occurrence.viewDependentCsgGeometry =
	viewDependentCsgGeometry ? TRUE : FALSE;
    occurrence.lodBacked = lodBacked ? TRUE : FALSE;
    occurrence.sourceMeshRequestValid = lodBacked ? TRUE : FALSE;
    if (lodBacked) {
	occurrence.sourceMeshRequest = sourceMeshRequest;
	compact_source_mesh_request_sync(occurrence.sourceMeshRequest,
	    occurrence.summary);
	compact_summary_lod_from_source_mesh_request(occurrence.summary,
	    occurrence.sourceMeshRequest);
    }
    occurrence.summary = compact_occurrence_summary(source,
	fullPath.c_str(), dp->d_namep,
	geometryKind && BU_STR_EQUAL(geometryKind, "annotation") ?
	"annotation" : typeLabel,
	geometryKind, revision, BObolRealizedShapeSummary::SHAPE_VLIST);
    occurrence.occurrenceIndex = source->occurrenceIndex.getValue();
    occurrence.booleanOperation = source->booleanOperation.getValue();
    occurrence.localTransform = pathMatrix;
    const int compacted = source->setCompactOccurrence(occurrence);
    if (!localBoundsValid && cadGeometry) {
	localBounds = compact_part_geometry_bounds(cadGeometry);
	localBoundsValid = !localBounds.isEmpty();
    }
    if (localBoundsValid)
	set_source_bounds_from_local_box(source, localBounds);
    if (compacted <= 0) {
	SbString msg;
	msg.sprintf("%s: direct compact wire geometry installation failed for primitive type '%s'",
	    fullPath.c_str(), typeLabel ? typeLabel : "");
	source->realizationDiagnostic = msg;
	return -1;
    }
    publish_direct_compact_occurrence(source, stream, dp, occurrence);
    return 1;
}

static SoBRLVListShape *
vlist_from_pnts(const struct rt_pnts_internal *pnts)
{
    if (!pnts || !pnts->point || pnts->count == 0)
	return NULL;
    RT_PNTS_CK_MAGIC(pnts);

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    std::vector<int> colorValid;
    std::vector<SbColor> colors;
    std::vector<int> scaleValid;
    std::vector<float> scales;
    std::vector<int> normalValid;
    std::vector<SbVec3f> normals;
    points.reserve(pnts->count);
    commands.reserve(pnts->count);

    colorValid.reserve(pnts->count);
    colors.reserve(pnts->count);
    scaleValid.reserve(pnts->count);
    scales.reserve(pnts->count);
    normalValid.reserve(pnts->count);
    normals.reserve(pnts->count);

    const double defaultScale = pnts->scale;
    auto appendPoint = [&](const fastf_t *v, const struct bu_color *c,
    const fastf_t *s, const fastf_t *n) {
	points.push_back(SbVec3f(static_cast<float>(v[X]),
				 static_cast<float>(v[Y]), static_cast<float>(v[Z])));
	commands.push_back(SoBRLVListShape::POINT);
	if (c) {
	    colorValid.push_back(1);
	    colors.push_back(SbColor(
				 static_cast<float>(c->buc_rgb[RED]),
				 static_cast<float>(c->buc_rgb[GRN]),
				 static_cast<float>(c->buc_rgb[BLU])));
	} else {
	    colorValid.push_back(0);
	    colors.push_back(SbColor(1.0f, 1.0f, 1.0f));
	}
	if (s && *s > 0.0) {
	    scaleValid.push_back(1);
	    scales.push_back(static_cast<float>(*s));
	} else if (!s && defaultScale > 0.0) {
	    scaleValid.push_back(1);
	    scales.push_back(static_cast<float>(defaultScale));
	} else {
	    scaleValid.push_back(0);
	    scales.push_back(0.0f);
	}
	if (n) {
	    normalValid.push_back(1);
	    normals.push_back(SbVec3f(static_cast<float>(n[X]),
				      static_cast<float>(n[Y]), static_cast<float>(n[Z])));
	} else {
	    normalValid.push_back(0);
	    normals.push_back(SbVec3f(0.0f, 0.0f, 1.0f));
	}
    };

    switch (pnts->type) {
	case RT_PNT_TYPE_PNT: {
	    const struct pnt *point = NULL;
	    for (BU_LIST_FOR(point, pnt, &(((struct pnt *)pnts->point)->l)))
		appendPoint(point->v, NULL, NULL, NULL);
	}
	break;
	case RT_PNT_TYPE_COL: {
	    const struct pnt_color *point = NULL;
	    for (BU_LIST_FOR(point, pnt_color, &(((struct pnt_color *)pnts->point)->l)))
		appendPoint(point->v, &point->c, NULL, NULL);
	}
	break;
	case RT_PNT_TYPE_SCA: {
	    const struct pnt_scale *point = NULL;
	    for (BU_LIST_FOR(point, pnt_scale, &(((struct pnt_scale *)pnts->point)->l)))
		appendPoint(point->v, NULL, &point->s, NULL);
	}
	break;
	case RT_PNT_TYPE_NRM: {
	    const struct pnt_normal *point = NULL;
	    for (BU_LIST_FOR(point, pnt_normal, &(((struct pnt_normal *)pnts->point)->l)))
		appendPoint(point->v, NULL, NULL, point->n);
	}
	break;
	case RT_PNT_TYPE_COL_SCA: {
	    const struct pnt_color_scale *point = NULL;
	    for (BU_LIST_FOR(point, pnt_color_scale, &(((struct pnt_color_scale *)pnts->point)->l)))
		appendPoint(point->v, &point->c, &point->s, NULL);
	}
	break;
	case RT_PNT_TYPE_COL_NRM: {
	    const struct pnt_color_normal *point = NULL;
	    for (BU_LIST_FOR(point, pnt_color_normal, &(((struct pnt_color_normal *)pnts->point)->l)))
		appendPoint(point->v, &point->c, NULL, point->n);
	}
	break;
	case RT_PNT_TYPE_SCA_NRM: {
	    const struct pnt_scale_normal *point = NULL;
	    for (BU_LIST_FOR(point, pnt_scale_normal, &(((struct pnt_scale_normal *)pnts->point)->l)))
		appendPoint(point->v, NULL, &point->s, point->n);
	}
	break;
	case RT_PNT_TYPE_COL_SCA_NRM: {
	    const struct pnt_color_scale_normal *point = NULL;
	    for (BU_LIST_FOR(point, pnt_color_scale_normal, &(((struct pnt_color_scale_normal *)pnts->point)->l)))
		appendPoint(point->v, &point->c, &point->s, point->n);
	}
	break;
	default:
	    return NULL;
    }

    if (points.empty() || points.size() != commands.size())
	return NULL;

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->setLineSet(points.data(), commands.data(), static_cast<int>(points.size()));
    shape->setPointAttributes(colorValid.data(), colors.data(),
			      scaleValid.data(), scales.data(),
			      normalValid.data(), normals.data(),
			      static_cast<int>(points.size()));
    return shape;
}

static int
cad_points_part_geometry_from_pnts(const struct rt_pnts_internal *pnts,
	Obol::PartGeometryBuilder &geometry)
{
    if (!pnts || !pnts->point || pnts->count == 0)
	return 0;
    RT_PNTS_CK_MAGIC(pnts);

    Obol::PointRep points;
    points.bounds.makeEmpty();
    points.positions.reserve(pnts->count);
    points.pointIds.reserve(pnts->count);
    points.colorValid.reserve(pnts->count);
    points.colors.reserve(pnts->count);
    points.scaleValid.reserve(pnts->count);
    points.scales.reserve(pnts->count);
    points.normalValid.reserve(pnts->count);
    points.normals.reserve(pnts->count);

    const double defaultScale = pnts->scale;
    auto appendPoint = [&](const fastf_t *v, const struct bu_color *c,
	const fastf_t *s, const fastf_t *n) {
	const SbVec3f point(static_cast<float>(v[X]),
	    static_cast<float>(v[Y]), static_cast<float>(v[Z]));
	points.positions.push_back(point);
	points.pointIds.push_back(
	    static_cast<uint32_t>(points.pointIds.size()));
	points.bounds.extendBy(point);
	if (c) {
	    points.colorValid.push_back(1u);
	    points.colors.push_back(SbColor(
		static_cast<float>(c->buc_rgb[RED]),
		static_cast<float>(c->buc_rgb[GRN]),
		static_cast<float>(c->buc_rgb[BLU])));
	} else {
	    points.colorValid.push_back(0u);
	    points.colors.push_back(SbColor(1.0f, 1.0f, 1.0f));
	}
	const float scale = s && *s > 0.0 ? static_cast<float>(*s) :
	    (!s && defaultScale > 0.0 ? static_cast<float>(defaultScale) :
	    0.0f);
	points.scaleValid.push_back(scale > 0.0f ? 1u : 0u);
	points.scales.push_back(scale);
	if (scale > 0.0f) {
	    const SbVec3f extent(scale, scale, scale);
	    points.bounds.extendBy(point - extent);
	    points.bounds.extendBy(point + extent);
	}
	if (n) {
	    points.normalValid.push_back(1u);
	    points.normals.push_back(SbVec3f(static_cast<float>(n[X]),
		static_cast<float>(n[Y]), static_cast<float>(n[Z])));
	} else {
	    points.normalValid.push_back(0u);
	    points.normals.push_back(SbVec3f(0.0f, 0.0f, 1.0f));
	}
    };

    switch (pnts->type) {
	case RT_PNT_TYPE_PNT: {
	    const struct pnt *point = NULL;
	    for (BU_LIST_FOR(point, pnt, &(((struct pnt *)pnts->point)->l)))
		appendPoint(point->v, NULL, NULL, NULL);
	}
	break;
	case RT_PNT_TYPE_COL: {
	    const struct pnt_color *point = NULL;
	    for (BU_LIST_FOR(point, pnt_color,
		&(((struct pnt_color *)pnts->point)->l)))
		appendPoint(point->v, &point->c, NULL, NULL);
	}
	break;
	case RT_PNT_TYPE_SCA: {
	    const struct pnt_scale *point = NULL;
	    for (BU_LIST_FOR(point, pnt_scale,
		&(((struct pnt_scale *)pnts->point)->l)))
		appendPoint(point->v, NULL, &point->s, NULL);
	}
	break;
	case RT_PNT_TYPE_NRM: {
	    const struct pnt_normal *point = NULL;
	    for (BU_LIST_FOR(point, pnt_normal,
		&(((struct pnt_normal *)pnts->point)->l)))
		appendPoint(point->v, NULL, NULL, point->n);
	}
	break;
	case RT_PNT_TYPE_COL_SCA: {
	    const struct pnt_color_scale *point = NULL;
	    for (BU_LIST_FOR(point, pnt_color_scale,
		&(((struct pnt_color_scale *)pnts->point)->l)))
		appendPoint(point->v, &point->c, &point->s, NULL);
	}
	break;
	case RT_PNT_TYPE_COL_NRM: {
	    const struct pnt_color_normal *point = NULL;
	    for (BU_LIST_FOR(point, pnt_color_normal,
		&(((struct pnt_color_normal *)pnts->point)->l)))
		appendPoint(point->v, &point->c, NULL, point->n);
	}
	break;
	case RT_PNT_TYPE_SCA_NRM: {
	    const struct pnt_scale_normal *point = NULL;
	    for (BU_LIST_FOR(point, pnt_scale_normal,
		&(((struct pnt_scale_normal *)pnts->point)->l)))
		appendPoint(point->v, NULL, &point->s, point->n);
	}
	break;
	case RT_PNT_TYPE_COL_SCA_NRM: {
	    const struct pnt_color_scale_normal *point = NULL;
	    for (BU_LIST_FOR(point, pnt_color_scale_normal,
		&(((struct pnt_color_scale_normal *)pnts->point)->l)))
		appendPoint(point->v, &point->c, &point->s, point->n);
	}
	break;
	default:
	    return 0;
    }

    if (points.positions.empty())
	return 0;
    geometry.points = std::move(points);
    return 1;
}


static SoBRLMeshShape *
mesh_from_bot(const struct rt_bot_internal *bot,
	      const SoBRLDatabaseSource *source)
{
    if (!bot || !bot->vertices || !bot->faces ||
	bot->num_vertices == 0 || bot->num_faces == 0 ||
	bot->num_vertices > INT_MAX || bot->num_faces > INT_MAX / 3)
	return NULL;
    RT_BOT_CK_MAGIC(bot);

    std::vector<SbVec3f> points;
    points.reserve(bot->num_vertices);
    for (size_t i = 0; i < bot->num_vertices; i++) {
	points.push_back(SbVec3f(static_cast<float>(bot->vertices[i * 3]),
				 static_cast<float>(bot->vertices[i * 3 + 1]),
				 static_cast<float>(bot->vertices[i * 3 + 2])));
    }

    std::vector<int32_t> indices;
    indices.reserve(bot->num_faces * 3);
    for (size_t i = 0; i < bot->num_faces; i++) {
	const int *face = &bot->faces[i * 3];
	if (face[0] < 0 || face[1] < 0 || face[2] < 0 ||
	    static_cast<size_t>(face[0]) >= bot->num_vertices ||
	    static_cast<size_t>(face[1]) >= bot->num_vertices ||
	    static_cast<size_t>(face[2]) >= bot->num_vertices)
	    return NULL;
	if (bot->orientation == RT_BOT_CW) {
	    indices.push_back(static_cast<int32_t>(face[0]));
	    indices.push_back(static_cast<int32_t>(face[2]));
	    indices.push_back(static_cast<int32_t>(face[1]));
	} else {
	    indices.push_back(static_cast<int32_t>(face[0]));
	    indices.push_back(static_cast<int32_t>(face[1]));
	    indices.push_back(static_cast<int32_t>(face[2]));
	}
    }

    std::vector<SbVec3f> normals;
    cad_bot_triangle_normals(normals, bot, points, indices);
    sanitize_triangle_normals(normals, points, indices);

    uint32_t threshold = source ? source->lodBotThreshold.getValue() : 0;
    SoBRLMeshShape *shape = (threshold > 0 &&
			     bot->num_faces >= static_cast<size_t>(threshold)) ?
			    new SoBRLLodMeshShape : new SoBRLMeshShape;
    shape->setIndexedTriangles(points.data(), static_cast<int>(points.size()),
			       indices.data(), static_cast<int>(indices.size()),
			       normals.empty() ? NULL : normals.data(),
			       static_cast<int>(normals.size()));
    return shape;
}


static void
primitive_indexed_face_set_free(struct rt_primitive_indexed_face_set *faceSet)
{
    rt_primitive_indexed_face_set_free(faceSet);
}

static int
indexed_face_finish(std::vector<int32_t> &face,
		    std::vector<SbVec3f> &faceNormals,
		    size_t pointCount,
		    std::vector<int32_t> &triangles,
		    std::vector<SbVec3f> *triangleNormals,
		    size_t *faceCount,
		    unsigned int *faceStamp,
		    std::vector<unsigned int> &seen)
{
    if (face.size() < 3)
	return 0;

    for (size_t i = 1; i + 1 < face.size(); i++) {
	triangles.push_back(face[0]);
	triangles.push_back(face[i]);
	triangles.push_back(face[i + 1]);
	if (triangleNormals && faceNormals.size() == face.size()) {
	    triangleNormals->push_back(faceNormals[0]);
	    triangleNormals->push_back(faceNormals[i]);
	    triangleNormals->push_back(faceNormals[i + 1]);
	}
    }

    face.clear();
    faceNormals.clear();
    if (faceCount)
	(*faceCount)++;
    if (faceStamp && seen.size() == pointCount) {
	if (*faceStamp == UINT_MAX) {
	    for (size_t i = 0; i < seen.size(); i++)
		seen[i] = 0;
	    *faceStamp = 1;
	} else {
	    (*faceStamp)++;
	}
    }
    return 1;
}

static int
indexed_faces_to_triangles(const int *indices,
			   size_t indexCount,
			   size_t pointCount,
			   std::vector<int32_t> &triangles,
			   const vect_t *normals = NULL,
			   size_t normalCount = 0,
			   std::vector<SbVec3f> *triangleNormals = NULL)
{
    if (!indices || !indexCount || !pointCount ||
	pointCount > static_cast<size_t>(INT_MAX) ||
	indexCount > static_cast<size_t>(INT_MAX))
	return 0;

    const int useNormals = normals && normalCount && triangleNormals;
    size_t normalIndex = 0;
    size_t faceCount = 0;
    unsigned int faceStamp = 1;
    std::vector<unsigned int> seen(pointCount, 0);
    std::vector<int32_t> face;
    std::vector<SbVec3f> faceNormals;

    for (size_t i = 0; i < indexCount; i++) {
	const int idx = indices[i];
	if (idx < 0) {
	    if (idx != -1 || !indexed_face_finish(face, faceNormals,
						  pointCount, triangles,
						  triangleNormals, &faceCount,
						  &faceStamp, seen))
		return 0;
	    continue;
	}

	if (static_cast<size_t>(idx) >= pointCount)
	    return 0;
	if (seen[static_cast<size_t>(idx)] == faceStamp)
	    return 0;
	seen[static_cast<size_t>(idx)] = faceStamp;
	face.push_back(static_cast<int32_t>(idx));
	if (useNormals) {
	    if (normalIndex >= normalCount)
		return 0;
	    faceNormals.push_back(SbVec3f(
				      static_cast<float>(normals[normalIndex][X]),
				      static_cast<float>(normals[normalIndex][Y]),
				      static_cast<float>(normals[normalIndex][Z])));
	    normalIndex++;
	}
    }

    if (!face.empty() && !indexed_face_finish(face, faceNormals, pointCount,
	    triangles, triangleNormals, &faceCount, &faceStamp, seen))
	return 0;
    if (useNormals && normalIndex != normalCount)
	return 0;
    if (useNormals && triangleNormals->size() != triangles.size())
	return 0;
    return faceCount > 0 && !triangles.empty();
}


static int
cad_mesh_part_geometry_from_indexed_face_set(
	const struct rt_primitive_indexed_face_set *faceSet,
	Obol::PartGeometryBuilder &geometry)
{
    if (!faceSet || !faceSet->points || !faceSet->point_count ||
	!faceSet->indices || !faceSet->index_count)
	return 0;

    size_t cornerCount = 0;
    for (size_t i = 0; i < faceSet->index_count; i++) {
	if (faceSet->indices[i] >= 0)
	    cornerCount++;
    }
    const int haveCompleteNormals =
	faceSet->normals && faceSet->normal_count == cornerCount;
    std::vector<int32_t> triangles;
    std::vector<SbVec3f> normals;
    if (!indexed_faces_to_triangles(faceSet->indices, faceSet->index_count,
	faceSet->point_count, triangles,
	haveCompleteNormals ? faceSet->normals : NULL,
	haveCompleteNormals ? faceSet->normal_count : 0,
	haveCompleteNormals ? &normals : NULL))
	return 0;

    Obol::TriMesh mesh;
    mesh.bounds.makeEmpty();
    mesh.positions.reserve(faceSet->point_count);
    for (size_t i = 0; i < faceSet->point_count; i++) {
	const SbVec3f point(
	    static_cast<float>(faceSet->points[i][X]),
	    static_cast<float>(faceSet->points[i][Y]),
	    static_cast<float>(faceSet->points[i][Z]));
	mesh.positions.push_back(point);
	mesh.bounds.extendBy(point);
    }
    /* Keep authored corner normals when present; otherwise publish smooth
     * crease-aware normals so Obol does not fall back to flat triangles. */
    sanitize_triangle_normals(normals, mesh.positions, triangles);
    mesh.indices.reserve(triangles.size());
    for (const int32_t index : triangles) {
	if (index < 0 || static_cast<size_t>(index) >= faceSet->point_count)
	    return 0;
	mesh.indices.push_back(static_cast<uint32_t>(index));
    }
    if (mesh.bounds.isEmpty() || mesh.indices.empty())
	return 0;
    if (!canonicalize_corner_normal_mesh(mesh, normals))
	return 0;
    geometry.shaded = std::move(mesh);
    return 1;
}

struct BObolOwnedStagedTriangleMesh {
    std::vector<fastf_t> points;
    std::vector<fastf_t> normals;
    std::vector<int> faces;
};

static int
cad_mesh_part_geometry_from_staged_source(
	const BObolStagedSourceMesh &staged,
	Obol::PartGeometryBuilder &geometry)
{
    if (!staged.isValid() ||
	staged.pointCount > static_cast<size_t>(UINT32_MAX) ||
	staged.faceCount > SIZE_MAX / 3)
	return 0;

    Obol::TriMesh mesh;
    mesh.bounds.makeEmpty();
    mesh.positions.reserve(staged.pointCount);
    for (size_t i = 0; i < staged.pointCount; ++i) {
	const SbVec3f point(
	    static_cast<float>(staged.points[i][X]),
	    static_cast<float>(staged.points[i][Y]),
	    static_cast<float>(staged.points[i][Z]));
	mesh.positions.push_back(point);
	mesh.bounds.extendBy(point);
    }

    const size_t indexCount = staged.faceCount * 3;
    std::vector<int32_t> indices;
    indices.reserve(indexCount);
    for (size_t i = 0; i < indexCount; ++i) {
	const int index = staged.faces[i];
	if (index < 0 || static_cast<size_t>(index) >= staged.pointCount)
	    return 0;
	indices.push_back(static_cast<int32_t>(index));
    }

    std::vector<SbVec3f> normals;
    if (staged.normals) {
	normals.reserve(indexCount);
	for (size_t i = 0; i < indexCount; ++i) {
	    normals.push_back(SbVec3f(
		static_cast<float>(staged.normals[i][X]),
		static_cast<float>(staged.normals[i][Y]),
		static_cast<float>(staged.normals[i][Z])));
	}
    }
    sanitize_triangle_normals(normals, mesh.positions, indices);

    mesh.indices.reserve(indexCount);
    for (const int32_t index : indices)
	mesh.indices.push_back(static_cast<uint32_t>(index));
    if (mesh.bounds.isEmpty() || mesh.indices.empty() ||
	!canonicalize_corner_normal_mesh(mesh, normals))
	return 0;

    geometry.shaded = std::move(mesh);
    geometry.shadedCullBackfaces = staged.shadedCullBackfaces != 0;
    return 1;
}

static unsigned long long
cad_brep_shaded_asset_key(struct db_i *dbip, struct directory *dp,
	const struct bg_tess_tol *ttol, const struct bn_tol *tol)
{
    if (!dbip || !dp || !ttol || !tol)
	return 0;
    struct bu_external external = BU_EXTERNAL_INIT_ZERO;
    if (db_get_external(&external, dp, dbip) != 0)
	return 0;
    struct bu_data_hash_state *hash = bu_data_hash_create();
    /* Version three requires a validated indexed-face-set whose referenced
     * triangles cover the source BREP boundary.  Earlier caches may contain
     * a successful subset after one or more CDT faces failed. */
    static const char contract[] =
	"BObol-display-asset:shaded-triangles:brep-indexed-face-set-v3:";
    bu_data_hash_update(hash, contract, sizeof(contract));
    bu_data_hash_update(hash, BOBOL_MESH_LOD_PROVIDER_VERSION,
	strlen(BOBOL_MESH_LOD_PROVIDER_VERSION));
    bu_data_hash_update(hash, external.ext_buf, external.ext_nbytes);
    const double tess[] = {
	ttol->abs, ttol->rel, ttol->norm, ttol->absmax, ttol->absmin,
	ttol->relmax, ttol->relmin, ttol->rel_lmax, ttol->rel_lmin
    };
    const double geometric[] = {
	tol->dist, tol->dist_sq, tol->perp, tol->para
    };
    bu_data_hash_update(hash, tess, sizeof(tess));
    bu_data_hash_update(hash, geometric, sizeof(geometric));
    unsigned long long key = bu_data_hash_val(hash);
    bu_data_hash_destroy(hash);
    bu_free_external(&external);
    return key ? key : 1;
}

static std::shared_ptr<BObolStagedSourceMesh>
cad_staged_mesh_from_primitive_face_set(
    struct db_i *dbip, struct directory *dp, struct rt_db_internal *intern,
    const struct bg_tess_tol *ttol, const struct bn_tol *tol,
    uint32_t sourceRevision, BObolSourceMeshRequest &request)
{
    request.clear();
    if (!dbip || !dp || !internal_payload_magic_valid(intern) ||
	!intern->idb_meth->ft_indexed_face_set || !ttol || !tol)
	return std::shared_ptr<BObolStagedSourceMesh>();

    struct rt_primitive_indexed_face_set faceSet;
    struct bv_view_info viewInfo = BV_VIEW_INFO_INIT;
    memset(&faceSet, 0, sizeof(faceSet));
    if (intern->idb_meth->ft_indexed_face_set(&faceSet, intern, ttol, tol,
	    &viewInfo) != BRLCAD_OK) {
	primitive_indexed_face_set_free(&faceSet);
	return std::shared_ptr<BObolStagedSourceMesh>();
    }

    size_t cornerCount = 0;
    for (size_t i = 0; i < faceSet.index_count; ++i)
	if (faceSet.indices[i] >= 0)
	    ++cornerCount;
    const bool haveNormals =
	faceSet.normals && faceSet.normal_count == cornerCount;
    std::vector<int32_t> triangles;
    std::vector<SbVec3f> triangleNormals;
    if (!indexed_faces_to_triangles(faceSet.indices, faceSet.index_count,
	    faceSet.point_count, triangles,
	    haveNormals ? faceSet.normals : NULL,
	    haveNormals ? faceSet.normal_count : 0,
	    haveNormals ? &triangleNormals : NULL) ||
	triangles.size() / 3 > static_cast<size_t>(INT_MAX)) {
	primitive_indexed_face_set_free(&faceSet);
	return std::shared_ptr<BObolStagedSourceMesh>();
    }

    std::shared_ptr<BObolOwnedStagedTriangleMesh> owned =
	std::make_shared<BObolOwnedStagedTriangleMesh>();
    owned->points.resize(faceSet.point_count * 3);
    SbBox3f bounds;
    bounds.makeEmpty();
    for (size_t i = 0; i < faceSet.point_count; ++i) {
	for (size_t axis = 0; axis < 3; ++axis) {
	    const fastf_t value = faceSet.points[i][axis];
	    if (!std::isfinite(value)) {
		primitive_indexed_face_set_free(&faceSet);
		return std::shared_ptr<BObolStagedSourceMesh>();
	    }
	    owned->points[i * 3 + axis] = value;
	}
	bounds.extendBy(SbVec3f(
	    static_cast<float>(faceSet.points[i][X]),
	    static_cast<float>(faceSet.points[i][Y]),
	    static_cast<float>(faceSet.points[i][Z])));
    }
    owned->faces.reserve(triangles.size());
    for (const int32_t index : triangles)
	owned->faces.push_back(static_cast<int>(index));
    if (haveNormals && triangleNormals.size() == triangles.size()) {
	owned->normals.resize(triangleNormals.size() * 3);
	for (size_t i = 0; i < triangleNormals.size(); ++i)
	    for (size_t axis = 0; axis < 3; ++axis)
		owned->normals[i * 3 + axis] = triangleNormals[i][axis];
    }
    primitive_indexed_face_set_free(&faceSet);
    if (bounds.isEmpty() || owned->faces.empty())
	return std::shared_ptr<BObolStagedSourceMesh>();

    std::shared_ptr<BObolStagedSourceMesh> staged =
	std::make_shared<BObolStagedSourceMesh>();
    staged->owner = owned;
    staged->points = reinterpret_cast<const point_t *>(owned->points.data());
    staged->normals = owned->normals.empty() ? NULL :
	reinterpret_cast<const vect_t *>(owned->normals.data());
    staged->faces = owned->faces.data();
    staged->pointCount = owned->points.size() / 3;
    staged->faceCount = owned->faces.size() / 3;
    staged->contentKey = cad_brep_shaded_asset_key(dbip, dp, ttol, tol);
    staged->shadedCullBackfaces = 0;
    staged->assetName = dp->d_namep ? dp->d_namep : "";
    staged->sourceRevision = sourceRevision;
    staged->byteCount = owned->points.size() * sizeof(fastf_t) +
	owned->normals.size() * sizeof(fastf_t) +
	owned->faces.size() * sizeof(int);
    if (!staged->contentKey || !staged->isValid())
	return std::shared_ptr<BObolStagedSourceMesh>();

    request.faceCount = staged->faceCount;
    request.pointCount = staged->pointCount;
    request.bounds = bounds;
    request.meshAssetBounds = bounds;
    request.meshAssetName = staged->assetName;
    request.meshAssetTessellationAbsTol = ttol->abs;
    request.meshAssetTessellationRelTol = ttol->rel;
    request.meshAssetTessellationNormTol = ttol->norm;
    return staged;
}

std::shared_ptr<BObolStagedSourceMesh>
bobol_database_brep_staged_mesh_variant(
    struct db_i *dbip, const char *name, const struct bg_tess_tol *ttol,
    uint64_t contentKey, uint32_t sourceRevision,
    BObolSourceMeshRequest &request)
{
    request.clear();
    if (!dbip || !name || !name[0] || !ttol || !contentKey)
	return std::shared_ptr<BObolStagedSourceMesh>();
    struct directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    if (!dp || dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BREP)
	return std::shared_ptr<BObolStagedSourceMesh>();
    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, dbip, NULL) < 0)
	return std::shared_ptr<BObolStagedSourceMesh>();
    const struct bn_tol tol = BN_TOL_INIT_TOL;
    std::shared_ptr<BObolStagedSourceMesh> staged =
	cad_staged_mesh_from_primitive_face_set(
	    dbip, dp, &intern, ttol, &tol, sourceRevision, request);
    rt_db_free_internal(&intern);
    if (!staged)
	return staged;
    /* The view allocator derives this identity from the validated canonical
     * asset plus a discrete tolerance band.  Do not let the helper's direct
     * database hash create a second logical identity for the same band. */
    staged->contentKey = contentKey;
    return staged;
}

static int
cad_mesh_part_geometry_from_primitive_face_set(
	struct rt_db_internal *intern, const struct bg_tess_tol *ttol,
	const struct bn_tol *tol, Obol::PartGeometryBuilder &geometry)
{
    if (!internal_payload_magic_valid(intern) ||
	!intern->idb_meth->ft_indexed_face_set || !ttol || !tol)
	return 0;

    struct rt_primitive_indexed_face_set faceSet;
    struct bv_view_info viewInfo = BV_VIEW_INFO_INIT;
    memset(&faceSet, 0, sizeof(faceSet));
    const int ret = intern->idb_meth->ft_indexed_face_set(&faceSet, intern,
	ttol, tol, &viewInfo);
    const int converted = ret == BRLCAD_OK ?
	cad_mesh_part_geometry_from_indexed_face_set(&faceSet, geometry) : 0;
    primitive_indexed_face_set_free(&faceSet);
    return converted;
}

static SoBRLMeshShape *
mesh_from_indexed_face_set(const struct rt_primitive_indexed_face_set *faceSet,
			   const SoBRLDatabaseSource *source)
{
    if (!faceSet || !faceSet->points || !faceSet->point_count ||
	!faceSet->indices || !faceSet->index_count ||
	faceSet->point_count > static_cast<size_t>(INT_MAX))
	return NULL;

    size_t cornerCount = 0;
    for (size_t i = 0; i < faceSet->index_count; i++) {
	if (faceSet->indices[i] >= 0)
	    cornerCount++;
    }
    const int haveCompleteNormals =
	faceSet->normals && faceSet->normal_count == cornerCount;
    std::vector<int32_t> triangles;
    std::vector<SbVec3f> normals;
    if (!indexed_faces_to_triangles(faceSet->indices, faceSet->index_count,
				    faceSet->point_count, triangles,
				    haveCompleteNormals ? faceSet->normals : NULL,
				    haveCompleteNormals ? faceSet->normal_count : 0,
				    haveCompleteNormals ? &normals : NULL))
	return NULL;
    if (triangles.size() > static_cast<size_t>(INT_MAX))
	return NULL;
    if (haveCompleteNormals && normals.size() > static_cast<size_t>(INT_MAX))
	return NULL;

    std::vector<SbVec3f> points;
    points.reserve(faceSet->point_count);
    for (size_t i = 0; i < faceSet->point_count; i++) {
	points.push_back(SbVec3f(
			     static_cast<float>(faceSet->points[i][X]),
			     static_cast<float>(faceSet->points[i][Y]),
			     static_cast<float>(faceSet->points[i][Z])));
    }
    sanitize_triangle_normals(normals, points, triangles);

    uint32_t threshold = source ? source->lodBotThreshold.getValue() : 0;
    SoBRLMeshShape *shape = (threshold > 0 &&
			     triangles.size() / 3 >= static_cast<size_t>(threshold)) ?
			    new SoBRLLodMeshShape : new SoBRLMeshShape;
    SbModernUtils::SoNodeRef owner(shape);
    shape->setIndexedTriangles(points.data(),
			       static_cast<int>(points.size()),
			       triangles.data(),
			       static_cast<int>(triangles.size()),
			       normals.empty() ? NULL : normals.data(),
			       static_cast<int>(normals.size()));
    owner.release()->unrefNoDelete();
    return shape;
}

static SoBRLMeshShape *
mesh_from_primitive_face_set(struct rt_db_internal *intern,
			     const SoBRLDatabaseSource *source)
{
    if (!internal_payload_magic_valid(intern) || !intern->idb_meth->ft_indexed_face_set)
	return NULL;

    struct rt_primitive_indexed_face_set faceSet;
    struct bv_view_info viewInfo = BV_VIEW_INFO_INIT;
    struct bg_tess_tol ttol = source_tess_tol(source);
    struct bn_tol tol = BN_TOL_INIT_TOL;
    memset(&faceSet, 0, sizeof(faceSet));

    int ret = intern->idb_meth->ft_indexed_face_set(&faceSet, intern,
	      &ttol, &tol, &viewInfo);
    if (ret != BRLCAD_OK) {
	primitive_indexed_face_set_free(&faceSet);
	return NULL;
    }

    SoBRLMeshShape *shape = mesh_from_indexed_face_set(&faceSet, source);
    primitive_indexed_face_set_free(&faceSet);
    return shape;
}

static void
publish_lod_result_metadata(SoBRLMeshShape *shape,
			    const BObolLodResult &result)
{
    if (!shape || result.providerStatus != BOBOL_LOD_PROVIDER_READY)
	return;

    (void)shape->applyStagedLodResult(result, &result.request);
}

static void
publish_lod_mesh_if_available(SoBRLMeshShape *shape,
			      const SoBRLDatabaseSource *source,
			      struct db_i *dbip,
			      const char *sourceName)
{
    if (!shape || !dbip || !sourceName)
	return;
    if (!source_mesh_lod_active(source))
	return;

    struct BObolMeshLodCacheStatus status =
	    BOBOL_MESH_LOD_CACHE_STATUS_INIT;
    if (bobol_mesh_lod_cache_status(dbip, sourceName, &status) != BRLCAD_OK)
	return;
    /* Scene realization is the foreground, time-to-first-frame path.  A
     * missing or stale LoD payload is generated by BObolLodService after the
     * source mesh request is published; generating it here makes a cold draw
     * block for tens of seconds on large BoTs. */
    if (!status.has_cache_key || !status.has_cached_payload ||
	status.stale_cache_entry)
	return;

    struct BObolMeshLod *lod = bobol_mesh_lod_get(dbip, sourceName);
    if (!lod)
	return;

    struct BObolMeshLodInfo info = BOBOL_MESH_LOD_INFO_INIT;
    struct BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    if (bobol_mesh_lod_hierarchy_info_get(lod, &hierarchy) &&
	bobol_mesh_lod_load_cut(lod, hierarchy.min_cut, 0) >= 0 &&
	bobol_mesh_lod_info_get(lod, &info)) {
	BObolLodRequest request;
	shape->makeLodRequest(request,
	    dbip->dbi_filename ? dbip->dbi_filename : "",
	    source ? source->sourceRevision.getValue() : 0,
	    source ? source->viewRevision.getValue() : 0,
	    0,
	    BOBOL_LOD_DRAW_SHADED,
	    "bobol_mesh_lod",
	    BOBOL_MESH_LOD_PROVIDER_VERSION,
	    BOBOL_LOD_QUALITY_FAST_DISPLAY);
	request.sourceContentHash = status.cache_key;
	if (request.objectPath.getLength() == 0)
	    request.objectPath = sourceName;
	if (request.objectName.getLength() == 0)
	    request.objectName = sourceName;
	request.bounds = SbBox3f(
			     SbVec3f(static_cast<float>(info.bmin[X]),
				     static_cast<float>(info.bmin[Y]),
				     static_cast<float>(info.bmin[Z])),
			     SbVec3f(static_cast<float>(info.bmax[X]),
				     static_cast<float>(info.bmax[Y]),
				     static_cast<float>(info.bmax[Z])));
	BObolLodResult result =
	    bobol_lod_result_from_mesh_lod_info(request, info, &status);
	struct BObolMeshLodData data;
	if (result.providerStatus == BOBOL_LOD_PROVIDER_READY &&
	    bobol_mesh_lod_data_get(lod, &data))
	    (void)bobol_lod_mesh_payload_from_mesh_lod_data(result.mesh,
		data);
	publish_lod_result_metadata(shape, result);
	bobol_mesh_lod_memshrink(lod);
    }

    bobol_mesh_lod_destroy(lod);
}

static SoBRLMeshShape *
mesh_from_nmg_region(struct nmgregion *r, struct bu_list *vlfree,
		     const struct bn_tol *tol)
{
    if (!r || !r->m_p || !vlfree || !tol)
	return NULL;

    /* librt's converter skips the quadratic edge-fusion pipeline for clean
     * primitive tessellations and retains it only as a degenerate fallback. */
    struct rt_bot_internal *bot = nmg_mdl_to_bot(r->m_p, vlfree, tol);
    if (!bot)
	return NULL;

    SoBRLMeshShape *shape = mesh_from_bot(bot, NULL);
    if (bot->vertices)
	bu_free(bot->vertices, "temporary tessellation BOT vertices");
    if (bot->faces)
	bu_free(bot->faces, "temporary tessellation BOT faces");
    bu_free(bot, "temporary tessellation BOT");
    return shape;
}

static SoBRLMeshShape *
mesh_from_tessellated_internal(struct rt_db_internal *intern,
			       const SoBRLDatabaseSource *source)
{
    if (!internal_payload_magic_valid(intern))
	return NULL;

    struct bg_tess_tol ttol = source_tess_tol(source);
    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct bu_list vlfree;
    struct model *m = NULL;
    struct nmgregion *r = NULL;
    SoBRLMeshShape *shape = NULL;
    int ret = -1;

    BU_LIST_INIT(&vlfree);
    m = nmg_mm();
    if (!m)
	return NULL;

    if (!BU_SETJUMP) {
	ret = rt_obj_tess(&r, m, intern, &ttol, &tol);
    } else {
	BU_UNSETJUMP;
	nmg_km(m);
	bu_list_free(&vlfree);
	return NULL;
    }
    BU_UNSETJUMP;

    if (ret == 0 && r)
	shape = mesh_from_nmg_region(r, &vlfree, &tol);

    nmg_km(m);
    bu_list_free(&vlfree);
    return shape;
}

static SoBRLMeshShape *
mesh_from_internal(struct rt_db_internal *intern,
		   const SoBRLDatabaseSource *source)
{
    if (!internal_payload_magic_valid(intern))
	return NULL;

    switch (intern->idb_type) {
	case ID_BOT:
	    return mesh_from_bot(static_cast<const struct rt_bot_internal *>(intern->idb_ptr), source);
	default:
	    break;
    }

    SoBRLMeshShape *faceSetShape = mesh_from_primitive_face_set(intern,
				   source);
    if (faceSetShape)
	return faceSetShape;

    return mesh_from_tessellated_internal(intern, source);
}

struct compact_mesh_prefill_job {
    compact_mesh_prefill_job(void) :
	ownsInternal(false),
	lodBacked(false),
	geometryTransform(SbMatrix::identity()),
	representative(NULL),
	success(false)
    {
	RT_DB_INTERNAL_INIT(&intern);
    }

    ~compact_mesh_prefill_job(void)
    {
	if (ownsInternal)
	    rt_db_free_internal(&intern);
    }

    struct directory *dp = NULL;
    struct rt_db_internal intern;
    bool ownsInternal;
    std::string cacheKey;
    std::string path;
    std::string sourceType;
    Obol::PartGeometryBuilder geometry;
    bool lodBacked;
    BObolSourceMeshRequest sourceMeshRequest;
    /* A non-null representative means this job reuses its immutable geometry
     * under geometryTransform rather than building a second mesh payload. */
    SbMatrix geometryTransform;
    compact_mesh_prefill_job *representative;
    bool success;
};

struct compact_mesh_prefill_collect {
    SoBRLDatabaseSource *source = NULL;
    BObolDatabaseSourceRealizationCache *cache = NULL;
    std::vector<std::unique_ptr<compact_mesh_prefill_job>> jobs;
    std::unordered_set<std::string> seen;
    SbString diagnostic;
    int failed = 0;
    /* Optional cancel handle: when set and cancelled, the walk stops gathering
     * new leaves so a cancelled deferred realization job aborts promptly. */
    const BObolCompactOccurrenceStream *cancel = NULL;
};

/* True when a streaming realization has been cancelled -- used to bail out of the
 * prefill's parallel phases promptly instead of finishing the whole batch. */
static inline bool
compact_mesh_prefill_cancelled(const BObolCompactOccurrenceStream *cancel)
{
    return cancel && cancel->isCancelled();
}

static union tree *
compact_mesh_prefill_collect_leaf(struct db_tree_state *tsp,
	const struct db_full_path *pathp, struct directory *dp,
	void *clientData)
{
    compact_mesh_prefill_collect *collect =
	static_cast<compact_mesh_prefill_collect *>(clientData);
    if (!collect || !collect->source || !collect->cache || !tsp ||
	!tsp->ts_dbip || !pathp || !dp)
	return TREE_NULL;

    if (compact_mesh_prefill_cancelled(collect->cancel))
	return make_nop_tree();

    SbBox3f cacheBounds;
    std::string cacheKey = realize_geometry_cache_key(dp);
    source_lod_cache_key_append(cacheKey, collect->source, cacheBounds,
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP);
    if (!collect->seen.insert(cacheKey).second)
	return make_nop_tree();
    if (collect->cache->findMeshVListCadGeometry(cacheKey) ||
	collect->cache->findMeshCadGeometry(cacheKey))
	return make_nop_tree();

    /* Only record the unique leaf here; the expensive rt_db_get_internal import
     * (and the type-dependent wire/special filter that needs it) is deferred to
     * a parallel pass in compact_mesh_prefill_cache, since concurrent reads of
     * the in-memory realization snapshot are lock-free. */
    std::unique_ptr<compact_mesh_prefill_job> job(
	new compact_mesh_prefill_job);
    job->dp = dp;
    job->cacheKey = cacheKey;
    char *path = db_path_to_string(pathp);
    if (path) {
	job->path = path;
	bu_free(path, "compact mesh prefill path");
    }
    collect->jobs.push_back(std::move(job));
    return make_nop_tree();
}

static const struct rt_bot_internal *
compact_mesh_prefill_bot(const compact_mesh_prefill_job &job)
{
    if (job.intern.idb_type != ID_BOT || !job.intern.idb_ptr)
	return NULL;
    return static_cast<const struct rt_bot_internal *>(job.intern.idb_ptr);
}

static bool
database_source_terminal_empty_bot(const struct rt_db_internal *intern)
{
    if (!intern || intern->idb_type != ID_BOT || !intern->idb_ptr)
	return false;
    const struct rt_bot_internal *bot =
	static_cast<const struct rt_bot_internal *>(intern->idb_ptr);
    /* A zero-face BoT has no drawable surface.  It is a valid no-op leaf,
     * commonly produced when a partial conversion cannot facetize one source
     * component.  A face-bearing BoT with missing vertices remains malformed
     * and follows the ordinary hard-failure path. */
    return bot->num_faces == 0;
}

static bool
compact_mesh_prefill_bot_cacheable(const struct rt_bot_internal *bot)
{
    /* Authored corner normals are geometry data.  A future extension can
     * transform and verify them, but never substitute them silently. */
    return bot && bot->vertices && bot->faces && bot->num_vertices > 0 &&
	bot->num_faces > 0 && !(bot->bot_flags & RT_BOT_HAS_SURFACE_NORMALS);
}

static bool
compact_mesh_prefill_bot_semantics_match(const struct rt_bot_internal *first,
	const struct rt_bot_internal *second)
{
    return compact_mesh_prefill_bot_cacheable(first) &&
	compact_mesh_prefill_bot_cacheable(second) &&
	first->mode == second->mode && first->orientation == second->orientation;
}


static bool
compact_mesh_prefill_pca_bucket(std::array<int64_t, 3> &bucket,
	const struct bg_pca_frame &frame)
{
    for (size_t i = 0; i < bucket.size(); i++) {
	const double value = static_cast<double>(frame.singular_values[i]);
	if (value <= SMALL_FASTF) {
	    bucket[i] = std::numeric_limits<int64_t>::min();
	    continue;
	}
	const double scaled = log2(value) * 4096.0;
	if (!std::isfinite(scaled) ||
	    scaled > static_cast<double>(std::numeric_limits<int64_t>::max()) ||
	    scaled < static_cast<double>(std::numeric_limits<int64_t>::min()))
	    return false;
	bucket[i] = static_cast<int64_t>(llround(scaled));
    }
    return true;
}

static std::string
compact_mesh_prefill_pca_bucket_key(const struct rt_bot_internal *bot,
	const std::array<int64_t, 3> &bucket)
{
    char key[256] = {0};
    snprintf(key, sizeof(key), "%zu:%zu:%u:%u:%" PRId64 ":%" PRId64 ":%" PRId64,
	bot->num_vertices, bot->num_faces, static_cast<unsigned int>(bot->mode),
	static_cast<unsigned int>(bot->orientation), bucket[0], bucket[1],
	bucket[2]);
    return key;
}

/* Per-job PCA data, precomputed in parallel (phase 1) so the serial match pass
 * (phase 2) never recomputes a signature. */
struct compact_mesh_prefill_pca {
    struct bg_trimesh_pca_signature signature;
    std::array<int64_t, 3> bucket;
    fastf_t comparisonTolerance;
    bool eligible;
};

static void
compact_mesh_prefill_find_transformed_reuse(
	std::vector<std::unique_ptr<compact_mesh_prefill_job>> &jobs,
	const BObolCompactOccurrenceStream *cancel)
{
    const size_t jobCount = jobs.size();
    if (!jobCount)
	return;

    /* Phase 1 (parallel): PCA signature + bucket per job.  Each entry is written
     * by exactly one thread and bg_trimesh_pca_get_signature only reads its
     * input mesh, so this parallelizes safely across distinct jobs. */
    std::vector<compact_mesh_prefill_pca> pca(jobCount);
    const auto computeSignature = [&](size_t i) {
	compact_mesh_prefill_pca &entry = pca[i];
	entry.eligible = false;
	entry.comparisonTolerance = VUNITIZE_TOL;
	const compact_mesh_prefill_job *job = jobs[i].get();
	const struct rt_bot_internal *bot = job ?
	    compact_mesh_prefill_bot(*job) : NULL;
	if (!compact_mesh_prefill_bot_cacheable(bot))
	    return;
	if (bg_trimesh_pca_get_signature(&entry.signature, bot->faces,
		bot->num_faces, reinterpret_cast<const point_t *>(bot->vertices),
		bot->num_vertices, VUNITIZE_TOL, 1.0e-6) != BRLCAD_OK)
	    return;
	if (!compact_mesh_prefill_pca_bucket(entry.bucket, entry.signature.frame))
	    return;
	/* PCA accumulation and a baked xpush transform both incur floating-point
	 * error proportional to object scale.  VUNITIZE_TOL alone is far too
	 * strict for million-unit meshes.  Validate every vertex and identical
	 * topology, but allow one part per billion of characteristic extent. */
	const double characteristicExtent =
	    static_cast<double>(entry.signature.frame.singular_values[0]) /
	    sqrt(static_cast<double>(bot->num_vertices));
	if (std::isfinite(characteristicExtent) && characteristicExtent > 0.0)
	    entry.comparisonTolerance = static_cast<fastf_t>(std::max(
		static_cast<double>(VUNITIZE_TOL),
		characteristicExtent * 1.0e-9));
	entry.eligible = true;
    };

    size_t threadCount = bu_avail_cpus();
    if (threadCount < 1)
	threadCount = 1;
    if (threadCount > jobCount)
	threadCount = jobCount;
    if (threadCount > 1) {
	std::atomic<size_t> cursor(0);
	auto poolWorker = [&]() {
	    for (size_t i = cursor.fetch_add(1); i < jobCount;
		 i = cursor.fetch_add(1)) {
		if (compact_mesh_prefill_cancelled(cancel))
		    return;
		computeSignature(i);
	    }
	};
	std::vector<std::thread> pool;
	pool.reserve(threadCount - 1);
	for (size_t t = 0; t + 1 < threadCount; t++)
	    pool.emplace_back(poolWorker);
	poolWorker();
	for (std::thread &thread : pool)
	    thread.join();
    } else {
	for (size_t i = 0; i < jobCount &&
	     !compact_mesh_prefill_cancelled(cancel); i++)
	    computeSignature(i);
    }
    if (compact_mesh_prefill_cancelled(cancel))
	return;

    /* Phase 2 (serial): match each eligible job against earlier representatives
     * using the precomputed signatures; identity of a representative is its job
     * index.  Processing in job order keeps the result deterministic and
     * identical to the former single-pass form. */
    std::unordered_map<unsigned long long, std::vector<size_t>> candidates;
    std::unordered_map<std::string, std::vector<size_t>> broadCandidates;
    for (size_t i = 0; i < jobCount; i++) {
	if (!pca[i].eligible)
	    continue;
	compact_mesh_prefill_job *job = jobs[i].get();
	const struct rt_bot_internal *bot = compact_mesh_prefill_bot(*job);
	const struct bg_trimesh_pca_signature &signature = pca[i].signature;
	const std::array<int64_t, 3> &bucket = pca[i].bucket;

	const auto matchesCandidate = [&](size_t candidateJob) {
	    compact_mesh_prefill_job *representative = jobs[candidateJob].get();
	    const struct rt_bot_internal *representativeBot =
		compact_mesh_prefill_bot(*representative);
	    const struct bg_trimesh_pca_signature &representativeSignature =
		pca[candidateJob].signature;
	    if (representative->lodBacked != job->lodBacked ||
		!compact_mesh_prefill_bot_semantics_match(representativeBot, bot) ||
		bg_trimesh_pca_equal(&representativeSignature,
		    representativeBot->faces, representativeBot->num_faces,
		    reinterpret_cast<const point_t *>(representativeBot->vertices),
		    representativeBot->num_vertices, &signature, bot->faces,
		    bot->num_faces, reinterpret_cast<const point_t *>(bot->vertices),
		    bot->num_vertices, std::max(pca[candidateJob].
			comparisonTolerance, pca[i].comparisonTolerance)) != 0)
		return false;

	    mat_t representativeToCandidate;
	    if (bg_pca_frame_relative_matrix(representativeToCandidate,
		&representativeSignature.frame, &signature.frame) != BRLCAD_OK)
		return false;
	    job->representative = representative;
	    job->geometryTransform = mat_to_sbmatrix(representativeToCandidate);
	    return true;
	};

	bool reused = false;
	const auto found = candidates.find(signature.hash);
	if (found != candidates.end()) {
	    for (size_t candidateJob : found->second) {
		if (matchesCandidate(candidateJob)) {
		    reused = true;
		    break;
		}
	    }
	}
	if (!reused) {
	    for (int xoffset = -1; xoffset <= 1 && !reused; xoffset++) {
		for (int yoffset = -1; yoffset <= 1 && !reused; yoffset++) {
		    for (int zoffset = -1; zoffset <= 1 && !reused; zoffset++) {
			std::array<int64_t, 3> nearby = bucket;
			const int offsets[3] = {xoffset, yoffset, zoffset};
			bool valid = true;
			for (size_t axis = 0; axis < nearby.size(); axis++) {
			    if (nearby[axis] == std::numeric_limits<int64_t>::min())
				continue;
			    if ((offsets[axis] > 0 && nearby[axis] >
				std::numeric_limits<int64_t>::max() - offsets[axis]) ||
				(offsets[axis] < 0 && nearby[axis] <
				std::numeric_limits<int64_t>::min() - offsets[axis])) {
				valid = false;
				break;
			    }
			    nearby[axis] += offsets[axis];
			}
			if (!valid)
			    continue;
			const auto broadFound = broadCandidates.find(
			    compact_mesh_prefill_pca_bucket_key(bot, nearby));
			if (broadFound == broadCandidates.end())
			    continue;
			for (size_t candidateJob : broadFound->second) {
			    if (matchesCandidate(candidateJob)) {
				reused = true;
				break;
			    }
			}
		    }
		}
	    }
	}
	if (reused)
	    continue;

	candidates[signature.hash].push_back(i);
	broadCandidates[compact_mesh_prefill_pca_bucket_key(bot, bucket)].push_back(
	    i);
    }
}

static void
compact_mesh_prefill_release_internal(compact_mesh_prefill_job &job)
{
    if (!job.ownsInternal)
	return;
    rt_db_free_internal(&job.intern);
    RT_DB_INTERNAL_INIT(&job.intern);
    job.ownsInternal = false;
}

/* Importing every xpush-expanded copy before matching defeats reuse on the
 * models where it matters most: the temporary internals alone may exceed
 * memory.  This bounded path imports one candidate at a time, retains only
 * unmatched static representatives, and releases every LoD-backed internal
 * once its canonical source identity and rigid placement are known. */
static void
compact_mesh_prefill_import_filter_reuse_bounded(
	compact_mesh_prefill_collect &collect, struct db_i *dbip)
{
    std::vector<std::unique_ptr<compact_mesh_prefill_job>> &jobs =
	collect.jobs;
    if (jobs.empty() || !dbip)
	return;

    const int drawMode = source_record_draw_mode(collect.source);
    std::vector<std::unique_ptr<compact_mesh_prefill_job>> kept;
    kept.reserve(jobs.size());
    std::vector<compact_mesh_prefill_pca> pca;
    pca.reserve(jobs.size());
    std::unordered_map<unsigned long long, std::vector<size_t>> candidates;
    std::unordered_map<std::string, std::vector<size_t>> broadCandidates;

    for (std::unique_ptr<compact_mesh_prefill_job> &jobPtr : jobs) {
	if (compact_mesh_prefill_cancelled(collect.cancel))
	    break;
	compact_mesh_prefill_job &imported = *jobPtr;
	if (rt_db_get_internal(&imported.intern, imported.dp, dbip, NULL) < 0 ||
	    !internal_payload_magic_valid(&imported.intern)) {
	    collect.failed = 1;
	    if (collect.diagnostic.getLength() == 0)
		collect.diagnostic.sprintf(
		    "%s: compact mesh internal fetch failed",
		    imported.path.empty() ?
			(imported.dp ? imported.dp->d_namep : "?") :
			imported.path.c_str());
	    continue;
	}
	imported.ownsInternal = true;
	imported.sourceType = primitive_type_label(&imported.intern);

	const int internalType = imported.intern.idb_type;
	if (database_source_terminal_empty_bot(&imported.intern)) {
	    compact_mesh_prefill_release_internal(imported);
	    continue;
	}
	const bool wireInstead =
	    (drawMode == BOBOL_LOD_DRAW_WIRE ||
	     (drawMode == BOBOL_LOD_DRAW_SHADED_BOTS &&
	      (!imported.intern.idb_meth ||
	       !imported.intern.idb_meth->ft_indexed_face_set))) &&
	    internalType != ID_BOT;
	const bool special = internalType == ID_MATERIAL ||
	    internalType == ID_PNTS ||
	    primitive_uses_wire_in_mesh_mode(internalType);
	if (wireInstead || special) {
	    compact_mesh_prefill_release_internal(imported);
	    continue;
	}
	if (internalType == ID_BOT &&
	    collect.source->lodBotThreshold.getValue() > 0) {
	    const struct rt_bot_internal *bot =
		static_cast<const struct rt_bot_internal *>(
		    imported.intern.idb_ptr);
	    if (bot && bot->num_faces >=
		collect.source->lodBotThreshold.getValue() &&
		cad_source_mesh_request_from_bot(imported.sourceMeshRequest, bot) &&
		cad_wire_part_geometry_from_aabb(
		    imported.sourceMeshRequest.bounds, imported.geometry)) {
		imported.sourceMeshRequest.meshAssetPath =
		    imported.path.c_str();
		imported.sourceMeshRequest.meshAssetName =
		    imported.dp && imported.dp->d_namep ?
			imported.dp->d_namep : "";
		imported.lodBacked = true;
		imported.success = true;
	    }
	}

	kept.push_back(std::move(jobPtr));
	const size_t i = kept.size() - 1;
	compact_mesh_prefill_job *job = kept[i].get();
	pca.emplace_back();
	compact_mesh_prefill_pca &entry = pca.back();
	entry.eligible = false;
	entry.comparisonTolerance = VUNITIZE_TOL;
	const struct rt_bot_internal *bot = compact_mesh_prefill_bot(*job);
	if (!compact_mesh_prefill_bot_cacheable(bot) ||
	    bg_trimesh_pca_get_signature(&entry.signature, bot->faces,
		bot->num_faces,
		reinterpret_cast<const point_t *>(bot->vertices),
		bot->num_vertices, VUNITIZE_TOL, 1.0e-6) != BRLCAD_OK ||
	    !compact_mesh_prefill_pca_bucket(entry.bucket,
		entry.signature.frame)) {
	    if (job->lodBacked)
		compact_mesh_prefill_release_internal(*job);
	    continue;
	}
	const double characteristicExtent =
	    static_cast<double>(entry.signature.frame.singular_values[0]) /
	    sqrt(static_cast<double>(bot->num_vertices));
	if (std::isfinite(characteristicExtent) && characteristicExtent > 0.0)
	    entry.comparisonTolerance = static_cast<fastf_t>(std::max(
		static_cast<double>(VUNITIZE_TOL),
		characteristicExtent * 1.0e-9));
	entry.eligible = true;

	const auto matchesCandidate = [&](size_t candidateJob) {
	    compact_mesh_prefill_job *representative =
		kept[candidateJob].get();
	    bool releaseRepresentative = false;
	    if (representative && !representative->ownsInternal) {
		if (rt_db_get_internal(&representative->intern,
			representative->dp, dbip, NULL) < 0 ||
		    !internal_payload_magic_valid(&representative->intern))
		    return false;
		representative->ownsInternal = true;
		releaseRepresentative = true;
	    }
	    const struct rt_bot_internal *representativeBot =
		compact_mesh_prefill_bot(*representative);
	    bool matched = representativeBot &&
		representative->lodBacked == job->lodBacked &&
		compact_mesh_prefill_bot_semantics_match(representativeBot,
		    bot) &&
		bg_trimesh_pca_equal(&pca[candidateJob].signature,
		    representativeBot->faces, representativeBot->num_faces,
		    reinterpret_cast<const point_t *>(
			representativeBot->vertices),
		    representativeBot->num_vertices, &entry.signature,
		    bot->faces, bot->num_faces,
		    reinterpret_cast<const point_t *>(bot->vertices),
		    bot->num_vertices, std::max(
			pca[candidateJob].comparisonTolerance,
			entry.comparisonTolerance)) == 0;
	    if (!matched) {
		if (releaseRepresentative)
		    compact_mesh_prefill_release_internal(*representative);
		return false;
	    }
	    mat_t representativeToCandidate;
	    matched = bg_pca_frame_relative_matrix(representativeToCandidate,
		&pca[candidateJob].signature.frame,
		&entry.signature.frame) == BRLCAD_OK;
	    if (releaseRepresentative)
		compact_mesh_prefill_release_internal(*representative);
	    if (!matched)
		return false;
	    job->representative = representative;
	    job->geometryTransform =
		mat_to_sbmatrix(representativeToCandidate);
	    return true;
	};

	bool reused = false;
	const auto found = candidates.find(entry.signature.hash);
	if (found != candidates.end()) {
	    for (size_t candidateJob : found->second) {
		if (matchesCandidate(candidateJob)) {
		    reused = true;
		    break;
		}
	    }
	}
	if (!reused) {
	    for (int xoffset = -1; xoffset <= 1 && !reused; xoffset++) {
		for (int yoffset = -1; yoffset <= 1 && !reused; yoffset++) {
		    for (int zoffset = -1; zoffset <= 1 && !reused;
			 zoffset++) {
			std::array<int64_t, 3> nearby = entry.bucket;
			const int offsets[3] = {xoffset, yoffset, zoffset};
			for (size_t axis = 0; axis < nearby.size(); axis++)
			    if (nearby[axis] !=
				std::numeric_limits<int64_t>::min())
				nearby[axis] += offsets[axis];
			const auto broadFound = broadCandidates.find(
			    compact_mesh_prefill_pca_bucket_key(bot, nearby));
			if (broadFound == broadCandidates.end())
			    continue;
			for (size_t candidateJob : broadFound->second) {
			    if (matchesCandidate(candidateJob)) {
				reused = true;
				break;
			    }
			}
		    }
		}
	    }
	}
	if (reused) {
	    job->success = true;
	    compact_mesh_prefill_release_internal(*job);
	    continue;
	}
	candidates[entry.signature.hash].push_back(i);
	broadCandidates[compact_mesh_prefill_pca_bucket_key(
	    bot, entry.bucket)].push_back(i);
	if (job->lodBacked)
	    compact_mesh_prefill_release_internal(*job);
    }

    jobs.swap(kept);
    /* Progressive jobs have already produced their box and canonical source
     * request.  The managed provider will import only the one representative
     * whose PoP cache is actually requested. */
    for (std::unique_ptr<compact_mesh_prefill_job> &job : jobs)
	if (job && job->lodBacked)
	    compact_mesh_prefill_release_internal(*job);
}

struct compact_mesh_prefill_workers {
    std::vector<std::unique_ptr<compact_mesh_prefill_job>> *jobs = NULL;
    struct bg_tess_tol ttol = BG_TESS_TOL_INIT_ZERO;
    struct bn_tol tol = BN_TOL_INIT_TOL;
    const BObolCompactOccurrenceStream *cancel = NULL;
};

static int
cad_mesh_part_geometry_from_tessellated_internal(
	struct rt_db_internal *intern, const struct bg_tess_tol *ttol,
	const struct bn_tol *tol, Obol::PartGeometryBuilder &geometry)
{
    if (!internal_payload_magic_valid(intern) || !ttol || !tol)
	return 0;

    struct bu_list vlfree;
    struct model *model = nmg_mm();
    struct nmgregion *region = NULL;
    int ret = -1;
    BU_LIST_INIT(&vlfree);
    if (!model)
	return 0;

    if (!BU_SETJUMP) {
	ret = rt_obj_tess(&region, model, intern, ttol, tol);
    } else {
	BU_UNSETJUMP;
	nmg_km(model);
	bu_list_free(&vlfree);
	return 0;
    }
    BU_UNSETJUMP;

    struct rt_bot_internal *bot = NULL;
    if (ret == 0 && region) {
	if (!BU_SETJUMP) {
	    bot = nmg_mdl_to_bot(region->m_p, &vlfree, tol);
	} else {
	    BU_UNSETJUMP;
	    bot = NULL;
	}
	BU_UNSETJUMP;
    }
    const int converted = bot ?
	cad_mesh_part_geometry_from_bot(bot, geometry) : 0;
    if (bot) {
	if (bot->vertices)
	    bu_free(bot->vertices, "prefill tessellation BOT vertices");
	if (bot->faces)
	    bu_free(bot->faces, "prefill tessellation BOT faces");
	bu_free(bot, "prefill tessellation BOT");
    }
    nmg_km(model);
    bu_list_free(&vlfree);
    return converted;
}


static int
cad_mesh_part_geometry_from_internal(struct rt_db_internal *intern,
	const SoBRLDatabaseSource *source, Obol::PartGeometryBuilder &geometry)
{
    if (!internal_payload_magic_valid(intern))
	return 0;
    if (intern->idb_type == ID_BOT)
	return cad_mesh_part_geometry_from_bot(
	    static_cast<const struct rt_bot_internal *>(intern->idb_ptr), geometry);

    struct bg_tess_tol ttol = source_tess_tol(source);
    struct bn_tol tol = BN_TOL_INIT_TOL;
    if (intern->idb_meth && intern->idb_meth->ft_indexed_face_set &&
	cad_mesh_part_geometry_from_primitive_face_set(intern, &ttol, &tol,
	    geometry))
	return 1;
    return cad_mesh_part_geometry_from_tessellated_internal(intern, &ttol,
	&tol, geometry);
}

/* Convert one collected prefill job's in-memory internal to compact geometry.
 * Only reads job.intern and writes job.geometry/job.success -- no shared state,
 * so BOT jobs run safely in parallel.  NMG tessellation (rt_obj_tess, used for
 * non-BOT non-indexed primitives) is NOT thread-safe (BU_SETJUMP/NMG globals),
 * so those jobs must stay on one thread (see compact_mesh_prefill_realize). */
static void
compact_mesh_prefill_realize_job(compact_mesh_prefill_job &job,
	const compact_mesh_prefill_workers *workers)
{
    if (job.lodBacked && job.success)
	return;
    if (job.representative) {
	job.success = true;
	return;
    }
    if (job.intern.idb_type == ID_BOT) {
	job.success = cad_mesh_part_geometry_from_bot(
	    static_cast<const struct rt_bot_internal *>(job.intern.idb_ptr),
	    job.geometry) != 0;
	return;
    }
    if (job.intern.idb_meth && job.intern.idb_meth->ft_indexed_face_set)
	job.success = cad_mesh_part_geometry_from_primitive_face_set(
	    &job.intern, &workers->ttol, &workers->tol, job.geometry) != 0;
    if (!job.success)
	job.success = cad_mesh_part_geometry_from_tessellated_internal(
	    &job.intern, &workers->ttol, &workers->tol, job.geometry) != 0;
}

static void
compact_mesh_prefill_realize(compact_mesh_prefill_workers *workers)
{
    if (!workers || !workers->jobs)
	return;
    std::vector<std::unique_ptr<compact_mesh_prefill_job>> &jobs =
	*workers->jobs;

    /* BOT-copy jobs are pure in-memory transforms and parallelize safely;
     * everything else (NMG tessellation, indexed-face-set) stays serial for
     * thread-safety.  Representatives carry no work. */
    std::vector<size_t> parallelJobs;
    std::vector<size_t> serialJobs;
    parallelJobs.reserve(jobs.size());
    for (size_t i = 0; i < jobs.size(); i++) {
	compact_mesh_prefill_job &job = *jobs[i];
	if (!job.representative && job.intern.idb_type == ID_BOT)
	    parallelJobs.push_back(i);
	else
	    serialJobs.push_back(i);
    }

    size_t threadCount = bu_avail_cpus();
    if (threadCount < 1)
	threadCount = 1;
    if (threadCount > parallelJobs.size())
	threadCount = parallelJobs.size();

    if (threadCount > 1) {
	std::atomic<size_t> cursor(0);
	auto poolWorker = [&]() {
	    for (size_t k = cursor.fetch_add(1); k < parallelJobs.size();
		 k = cursor.fetch_add(1)) {
		if (compact_mesh_prefill_cancelled(workers->cancel))
		    return;
		compact_mesh_prefill_realize_job(*jobs[parallelJobs[k]],
		    workers);
	    }
	};
	std::vector<std::thread> pool;
	pool.reserve(threadCount - 1);
	for (size_t t = 0; t + 1 < threadCount; t++)
	    pool.emplace_back(poolWorker);
	poolWorker();
	for (std::thread &thread : pool)
	    thread.join();
    } else {
	for (size_t idx : parallelJobs) {
	    if (compact_mesh_prefill_cancelled(workers->cancel))
		return;
	    compact_mesh_prefill_realize_job(*jobs[idx], workers);
	}
    }

    for (size_t idx : serialJobs) {
	if (compact_mesh_prefill_cancelled(workers->cancel))
	    return;
	compact_mesh_prefill_realize_job(*jobs[idx], workers);
    }
}

/* Import each gathered unique leaf's internal in parallel (concurrent reads of
 * the in-memory realization snapshot are lock-free), then serially drop the
 * jobs the main walk handles instead (wire/special primitives) and flag any
 * import failure.  Imports are the bulk of the former serial collect. */
static void
compact_mesh_prefill_import_and_filter(compact_mesh_prefill_collect &collect,
	struct db_i *dbip)
{
    std::vector<std::unique_ptr<compact_mesh_prefill_job>> &jobs = collect.jobs;
    const size_t jobCount = jobs.size();
    if (!jobCount || !dbip)
	return;

    const auto importJob = [&](size_t i) {
	compact_mesh_prefill_job &job = *jobs[i];
	if (rt_db_get_internal(&job.intern, job.dp, dbip, NULL) >= 0 &&
	    internal_payload_magic_valid(&job.intern)) {
	    job.ownsInternal = true;
	    job.sourceType = primitive_type_label(&job.intern);
	}
    };

    size_t threadCount = bu_avail_cpus();
    if (threadCount < 1)
	threadCount = 1;
    if (threadCount > jobCount)
	threadCount = jobCount;
    if (threadCount > 1) {
	std::atomic<size_t> cursor(0);
	auto poolWorker = [&]() {
	    for (size_t i = cursor.fetch_add(1); i < jobCount;
		 i = cursor.fetch_add(1)) {
		if (compact_mesh_prefill_cancelled(collect.cancel))
		    return;
		importJob(i);
	    }
	};
	std::vector<std::thread> pool;
	pool.reserve(threadCount - 1);
	for (size_t t = 0; t + 1 < threadCount; t++)
	    pool.emplace_back(poolWorker);
	poolWorker();
	for (std::thread &thread : pool)
	    thread.join();
    } else {
	for (size_t i = 0; i < jobCount; i++) {
	    if (compact_mesh_prefill_cancelled(collect.cancel))
		return;
	    importJob(i);
	}
    }
    if (compact_mesh_prefill_cancelled(collect.cancel))
	return;

    /* A failed import is a hard failure (matches the former inline behavior).
     * Otherwise drop wire/special leaves -- the main walk realizes those as
     * wireframe/points so the prefill must not cache them as meshes. */
    const int drawMode = source_record_draw_mode(collect.source);
    std::vector<std::unique_ptr<compact_mesh_prefill_job>> kept;
    kept.reserve(jobCount);
    for (std::unique_ptr<compact_mesh_prefill_job> &jobPtr : jobs) {
	compact_mesh_prefill_job &job = *jobPtr;
	if (!job.ownsInternal) {
	    collect.failed = 1;
	    if (collect.diagnostic.getLength() == 0)
		collect.diagnostic.sprintf(
		    "%s: compact mesh internal fetch failed",
		    job.path.empty() ? (job.dp ? job.dp->d_namep : "?") :
		    job.path.c_str());
	    continue;
	}
	const int internalType = job.intern.idb_type;
	if (database_source_terminal_empty_bot(&job.intern))
	    continue;
	const bool wireInstead =
	    (drawMode == BOBOL_LOD_DRAW_WIRE ||
	     (drawMode == BOBOL_LOD_DRAW_SHADED_BOTS &&
	      (!job.intern.idb_meth ||
	       !job.intern.idb_meth->ft_indexed_face_set))) &&
	    internalType != ID_BOT;
	const bool special = internalType == ID_MATERIAL ||
	    internalType == ID_PNTS ||
	    primitive_uses_wire_in_mesh_mode(internalType);
	if (wireInstead || special)
	    continue;
	if (internalType == ID_BOT &&
	    collect.source->lodBotThreshold.getValue() > 0) {
	    const struct rt_bot_internal *bot =
		static_cast<const struct rt_bot_internal *>(job.intern.idb_ptr);
	    if (bot && bot->num_faces >=
		collect.source->lodBotThreshold.getValue() &&
		cad_source_mesh_request_from_bot(job.sourceMeshRequest, bot) &&
		cad_wire_part_geometry_from_aabb(job.sourceMeshRequest.bounds,
		    job.geometry)) {
		job.sourceMeshRequest.meshAssetPath = job.path.c_str();
		job.sourceMeshRequest.meshAssetName =
		    job.dp && job.dp->d_namep ? job.dp->d_namep : "";
		/* The compact realization keeps only a box until the managed
		 * service loads the view-selected PoP cut. */
		job.lodBacked = true;
		job.success = true;
	    }
	}
	kept.push_back(std::move(jobPtr));
    }
    jobs.swap(kept);
}

static int
compact_mesh_prefill_cache(SoBRLDatabaseSource *source,
	BObolDatabaseSourceRealizationCache *cache, const char *treeName,
	const BObolCompactOccurrenceStream *cancel)
{
    struct db_i *dbip = source ? source->getDatabase() : NULL;
    if (!source || !cache || !dbip || !treeName || !treeName[0])
	return 0;

    const int prefillTiming = getenv("BOBOL_DRAW_TIMING") ? 1 : 0;
    const int64_t collectStart = prefillTiming ? bu_gettime() : 0;
    compact_mesh_prefill_collect collect;
    collect.source = source;
    collect.cache = cache;
    collect.cancel = cancel;
    struct db_tree_state initState;
    db_init_db_tree_state(&initState, dbip);
    initState.ts_stop_at_regions = 0;
    const char *av[1] = {treeName};
    const int walkRet = db_walk_tree_leaf_instances(dbip, 1, av, 1,
	&initState, NULL, NULL, compact_mesh_prefill_collect_leaf, &collect);
    db_free_db_tree_state(&initState);
    if (walkRet < 0) {
	SbString diagnostic;
	diagnostic.sprintf("%s: compact mesh occurrence discovery failed",
	    treeName);
	source->realizationDiagnostic = diagnostic;
	return -1;
    }
    if (compact_mesh_prefill_cancelled(cancel))
	return -1;
    if (collect.jobs.empty())
	return 0;

    /* Parallel import of the gathered unique leaves (formerly serial in the
     * walk), then drop the wire/special leaves the main walk handles. */
    const int64_t importStart = prefillTiming ? bu_gettime() : 0;
    uint64_t encodedLeafBytes = 0;
    for (const std::unique_ptr<compact_mesh_prefill_job> &job : collect.jobs) {
	const uint64_t leafBytes = job && job->dp && job->dp->d_len > 0 ?
	    static_cast<uint64_t>(job->dp->d_len) : 0;
	encodedLeafBytes = UINT64_MAX - encodedLeafBytes < leafBytes ?
	    UINT64_MAX : encodedLeafBytes + leafBytes;
    }
    const bool boundedReuse =
	collect.jobs.size() > 1 &&
	encodedLeafBytes > 512ULL * 1024ULL * 1024ULL;
    if (boundedReuse)
	compact_mesh_prefill_import_filter_reuse_bounded(collect, dbip);
    else
	compact_mesh_prefill_import_and_filter(collect, dbip);
    if (compact_mesh_prefill_cancelled(cancel))
	return -1;
    if (collect.failed) {
	source->realizationDiagnostic = collect.diagnostic.getLength() > 0 ?
	    collect.diagnostic : SbString("compact mesh internal fetch failed");
	return -1;
    }
    if (collect.jobs.empty())
	return 0;

    const int64_t reuseStart = prefillTiming ? bu_gettime() : 0;
    if (!boundedReuse)
	compact_mesh_prefill_find_transformed_reuse(collect.jobs, cancel);
    if (compact_mesh_prefill_cancelled(cancel))
	return -1;
    if (prefillTiming) {
	size_t reused = 0;
	for (const std::unique_ptr<compact_mesh_prefill_job> &jobPtr :
	     collect.jobs)
	    if (jobPtr && jobPtr->representative)
		reused++;
	bu_log("[obol-timing] prefill: walk %.1f ms, import %.1f ms, "
	    "find_reuse %.1f ms; %zu unique jobs, %zu reused via PCA%s\n",
	    (double)(importStart - collectStart) / 1000.0,
	    (double)(reuseStart - importStart) / 1000.0,
	    (double)(bu_gettime() - reuseStart) / 1000.0,
	    collect.jobs.size(), reused,
	    boundedReuse ? " (bounded import)" : "");
    }
    compact_mesh_prefill_workers workers;
    workers.jobs = &collect.jobs;
    workers.ttol = source_tess_tol(source);
    workers.cancel = cancel;
    const int64_t realizeStart = prefillTiming ? bu_gettime() : 0;
    compact_mesh_prefill_realize(&workers);
    if (compact_mesh_prefill_cancelled(cancel))
	return -1;
    if (prefillTiming)
	bu_log("[obol-timing] prefill: realize %.1f ms\n",
	    (double)(bu_gettime() - realizeStart) / 1000.0);

    for (std::unique_ptr<compact_mesh_prefill_job> &job : collect.jobs) {
	if (!job->success) {
	    SbString diagnostic;
	    diagnostic.sprintf(
		"%s: compact mesh conversion/tessellation failed for primitive type '%s'",
		job->path.empty() ? job->dp->d_namep : job->path.c_str(),
		job->sourceType.c_str());
	    source->realizationDiagnostic = diagnostic;
	    return -1;
	}
	if (job->representative) {
	    const BObolCachedPartGeometry *representative =
		cache->findMeshCadGeometry(job->representative->cacheKey);
	    if (!representative || !representative->geometry) {
		source->realizationDiagnostic =
		    "compact mesh transformed-geometry representative is unavailable";
		return -1;
	    }
	    const SbBox3f bounds = database_source_transform_bounds(
		representative->bounds, job->geometryTransform);
	    const bool lodBacked = job->lodBacked &&
		representative->sourceMeshRequestValid;
	    BObolSourceMeshRequest sourceMeshRequest =
		representative->sourceMeshRequest;
	    sourceMeshRequest.meshAssetTransform = job->geometryTransform;
	    cache->storeMeshCadGeometryReference(job->cacheKey,
		representative->geometry, job->geometryTransform,
		job->sourceType.c_str(), lodBacked ? "aabb" : "surface",
		&bounds, lodBacked,
		lodBacked ? &sourceMeshRequest : NULL);
	    continue;
	}
	BObolSourceMeshRequest sourceMeshRequest = job->sourceMeshRequest;
	bool lodBacked = job->lodBacked;
	SbBox3f bounds;
	bounds.makeEmpty();
	if (job->geometry.points)
	    bounds.extendBy(job->geometry.points->bounds);
	if (job->geometry.wire)
	    bounds.extendBy(job->geometry.wire->bounds);
	if (job->geometry.shaded)
	    bounds.extendBy(job->geometry.shaded->bounds);
	cache->storeMeshCadGeometry(job->cacheKey, std::move(job->geometry),
	    job->sourceType.c_str(), lodBacked ? "aabb" : "surface",
	    &bounds, lodBacked,
	    lodBacked ? &sourceMeshRequest : NULL);
    }
    return static_cast<int>(collect.jobs.size());
}

static std::string
database_source_evaluated_path_string(const SoBRLDatabaseSource *source)
{
    if (!source)
	return std::string();
    return std::string(database_source_skip_leading_slash(
			   source->path.getValue().getString()));
}

static SoBRLVListShape *
vlist_from_evaluated_wire_path(SoBRLDatabaseSource *source)
{
    struct db_i *dbip = source ? source->getDatabase() : NULL;
    if (!source || !dbip)
	return NULL;

    const std::string path = database_source_evaluated_path_string(source);
    if (path.empty())
	return NULL;

    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct bg_tess_tol ttol = source_tess_tol(source);
    struct bu_list vhead;
    struct bu_list vlfree;
    BU_LIST_INIT(&vhead);
    BU_LIST_INIT(&vlfree);
    const auto releasePlot = [&](bu_list *head) {
	RT_FREE_VLIST(&vlfree, head);
	bg_vlist_cleanup(&vlfree);
    };
    std::unique_ptr<bu_list, decltype(releasePlot)> plot(&vhead, releasePlot);

    struct rt_eval_wireframe_opts opts = RT_EVAL_WIREFRAME_OPTS_INIT;
    const int ret = rt_eval_wireframe(&vhead, &vlfree, dbip, path.c_str(),
				      &tol, &ttol, &opts);
    if (ret != BRLCAD_OK)
	return NULL;

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    convert_vlist(points, commands, &vhead);
    if (points.empty() || points.size() != commands.size() ||
	points.size() > static_cast<size_t>(INT_MAX))
	return NULL;

    SoBRLVListShape *shape = new SoBRLVListShape;
    SbModernUtils::SoNodeRef owner(shape);
    shape->setLineSet(points.data(), commands.data(),
		      static_cast<int>(points.size()));
    owner.release()->unrefNoDelete();
    return shape;
}

static SoBRLMeshShape *
mesh_from_evaluated_points_path(SoBRLDatabaseSource *source,
	SbBox3f &sourceBounds)
{
    sourceBounds.makeEmpty();
    struct db_i *dbip = source ? source->getDatabase() : NULL;
    if (!source || !dbip)
	return NULL;

    const std::string path = database_source_evaluated_path_string(source);
    if (path.empty())
	return NULL;

    struct rt_primitive_indexed_face_set faceSet;
    memset(&faceSet, 0, sizeof(faceSet));
    std::unique_ptr<rt_primitive_indexed_face_set, decltype(&bobol_evaluated_points_face_set_free)>
	ownedFaceSet(&faceSet, bobol_evaluated_points_face_set_free);
    const int ret = bobol_evaluated_points_evaluate_path_face_set(
			dbip, path.c_str(), &faceSet);
    if (ret != BRLCAD_OK)
	return NULL;

    if (faceSet.source_bounds_valid) {
	sourceBounds = database_source_box_from_minmax(
	    SbVec3f(static_cast<float>(faceSet.source_bounds_min[X]),
		    static_cast<float>(faceSet.source_bounds_min[Y]),
		    static_cast<float>(faceSet.source_bounds_min[Z])),
	    SbVec3f(static_cast<float>(faceSet.source_bounds_max[X]),
		    static_cast<float>(faceSet.source_bounds_max[Y]),
		    static_cast<float>(faceSet.source_bounds_max[Z])));
    }

    return mesh_from_indexed_face_set(&faceSet, source);
}

static std::string
database_source_evaluated_display_name(const SoBRLDatabaseSource *source)
{
    if (!source)
	return std::string();

    std::string name = database_source_leaf_component(source->path.getValue());
    if (!name.empty())
	return name;
    return stable_name_from_path(source->path.getValue().getString(), 1);
}

static int
realize_evaluated_wire_source(SoBRLDatabaseSource *source, uint32_t revision)
{
    if (!source)
	return -1;

    SoBRLVListShape *shape = vlist_from_evaluated_wire_path(source);
    if (!shape) {
	source->realizationDiagnostic =
	    "evaluated-wire provider produced no drawable geometry";
	return -1;
    }

    SbModernUtils::SoNodeRef owner(shape);
    const std::string fullPath =
	database_source_full_path_string(source->path.getValue());
    const std::string displayName =
	database_source_evaluated_display_name(source);
    assign_realized_identity(shape, NULL, fullPath.c_str(),
			     displayName.c_str(), "evaluated-wire",
			     revision, source);
    shape->geometryKind = "evaluated-wire";
    source->setRealizationRoleFlags(SoBRLDatabaseSource::REALIZATION_ROLE_CSG);
    database_source_add_realized_child(source, shape);
    return 1;
}

static int
realize_evaluated_points_source(SoBRLDatabaseSource *source, uint32_t revision)
{
    if (!source)
	return -1;

    SbBox3f sourceBounds;
    SoBRLMeshShape *shape = mesh_from_evaluated_points_path(source,
						    sourceBounds);
    if (!shape) {
	source->realizationDiagnostic =
	    "evaluated-points provider produced no drawable geometry";
	return -1;
    }

    SbModernUtils::SoNodeRef owner(shape);
    const std::string fullPath =
	database_source_full_path_string(source->path.getValue());
    const std::string displayName =
	database_source_evaluated_display_name(source);
    assign_realized_identity(shape, NULL, fullPath.c_str(),
			     displayName.c_str(), "evaluated-points",
			     revision, source);
    shape->geometryKind = "evaluated-points";
    source->setRealizationRoleFlags(
	SoBRLDatabaseSource::REALIZATION_ROLE_CSG |
	SoBRLDatabaseSource::REALIZATION_ROLE_MESH);
    database_source_add_realized_child(source, shape);
    if (!sourceBounds.isEmpty())
	set_source_bounds_from_local_box(source, sourceBounds, TRUE);
    return 1;
}

static SoBRLMeshShape *
mesh_instance_for_shared_geometry(const SoBRLMeshShape *sharedShape)
{
    if (sharedShape &&
	sharedShape->isOfType(SoBRLLodMeshShape::getClassTypeId()))
	return new SoBRLLodMeshShape;
    return new SoBRLMeshShape;
}

static int
realize_direct_leaf_mesh(SoBRLDatabaseSource *source,
			 BObolDatabaseSourceRealizationCache *cache,
			 uint32_t revision)
{
    struct db_i *dbip = source ? source->getDatabase() : NULL;
    if (!source || !cache || !dbip)
	return 0;

    const std::string leafName =
	database_source_leaf_component(source->path.getValue());
    if (leafName.empty())
	return 0;

    struct directory *dp =
	db_lookup(dbip, leafName.c_str(), LOOKUP_QUIET);
    if (!dp || (dp->d_flags & RT_DIR_COMB))
	return 0;

    const std::string fullPath =
	database_source_full_path_string(source->path.getValue());
    SbMatrix pathMatrix;
    if (!database_source_path_matrix(dbip, source->path.getValue(),
	pathMatrix))
	return 0;
    SoBRLVListShape *sharedVListShape = NULL;
    SoBRLMeshShape *sharedMeshShape = NULL;
    SbBox3f cacheBounds;
    std::string cacheKey = realize_geometry_cache_key(dp);
    source_lod_cache_key_append(cacheKey, source, cacheBounds,
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP);
    BObolRealizationCacheMap<SoBRLVListShape *>::iterator foundVList =
	cache->sharedMeshVListGeometry.find(cacheKey);
    if (foundVList != cache->sharedMeshVListGeometry.end()) {
	sharedVListShape = foundVList->second;
    } else {
	BObolRealizationCacheMap<SoBRLMeshShape *>::iterator foundMesh =
	    cache->sharedMeshGeometry.find(cacheKey);
	if (foundMesh != cache->sharedMeshGeometry.end())
	    sharedMeshShape = foundMesh->second;
    }
    bobol_performance_counter_add(
	(sharedVListShape || sharedMeshShape) ? BOBOL_PERF_MESH_CACHE_HITS :
	BOBOL_PERF_MESH_CACHE_MISSES, 1);

    const char *typeLabel = sharedVListShape ?
			    sharedVListShape->sourceType.getValue().getString() :
			    (sharedMeshShape ?
			     sharedMeshShape->sourceType.getValue().getString() :
			     NULL);
    if (!sharedVListShape && !sharedMeshShape) {
	owned_leaf_internal validInternal;
	if (rt_db_get_internal(&validInternal.local, dp, dbip, NULL) < 0 ||
	    !internal_payload_magic_valid(&validInternal.local)) {
	    SbString msg;
	    msg.sprintf("%s: direct leaf mesh internal fetch failed",
			fullPath.c_str());
	    source->realizationDiagnostic = msg;
	    if (validInternal.local.idb_ptr)
		rt_db_free_internal(&validInternal.local);
	    return -1;
	}
	validInternal.ownsLocal = true;

	typeLabel = primitive_type_label(&validInternal.local);
	const int internalType = validInternal.local.idb_type;
	SbBox3f localBounds;
	if (source_view_lod_active(source))
	    (void)local_bounds_from_internal(&validInternal.local, localBounds);
	if (internalType == ID_MATERIAL) {
	    SoBRLMaterialObject *materialObject =
		material_object_from_internal(
		    static_cast<struct rt_material_internal *>(
			validInternal.local.idb_ptr));
	    if (!materialObject) {
		SbString msg;
		msg.sprintf("%s: material object realization failed",
			    fullPath.c_str());
		source->realizationDiagnostic = msg;
		return -1;
	    }
	    SoSeparator *leaf = realize_matrix_leaf_separator(pathMatrix);
	    assign_material_identity(materialObject,
				     fullPath.c_str(),
				     dp->d_namep, typeLabel, revision);
	    leaf->addChild(materialObject);
	    database_source_add_realized_child(source, leaf);
	    return 1;
	}

	if (internalType == ID_PNTS) {
	    sharedVListShape = vlist_from_pnts(
		static_cast<const struct rt_pnts_internal *>(
		    validInternal.local.idb_ptr));
	    if (sharedVListShape)
		assign_shared_geometry_identity(sharedVListShape,
		    dp->d_namep, typeLabel, revision, "point");
	} else if (primitive_uses_wire_in_mesh_mode(internalType)) {
	    sharedVListShape = vlist_from_lod_realization_internal(
		&validInternal.local, source, localBounds);
	    if (!sharedVListShape)
		sharedVListShape = vlist_from_plot_internal(&validInternal.local,
		    source);
	    if (sharedVListShape) {
		assign_shared_geometry_identity(sharedVListShape,
		    dp->d_namep, typeLabel, revision, "line");
		if (primitive_is_annotation(internalType, typeLabel)) {
		    sharedVListShape->sourceType = "annotation";
		    sharedVListShape->geometryKind = "annotation";
		}
	    }
	} else if ((source_record_draw_mode(source) == BOBOL_LOD_DRAW_WIRE ||
		    (source_record_draw_mode(source) ==
			BOBOL_LOD_DRAW_SHADED_BOTS &&
		     (!validInternal.local.idb_meth ||
		      !validInternal.local.idb_meth->ft_indexed_face_set))) &&
			   internalType != ID_BOT) {
	    sharedVListShape = vlist_from_lod_realization_internal(
		&validInternal.local, source, localBounds);
	    if (!sharedVListShape)
		sharedVListShape = vlist_from_plot_internal(&validInternal.local,
		    source);
	    if (sharedVListShape)
		assign_shared_geometry_identity(sharedVListShape,
		    dp->d_namep, typeLabel, revision, "line");
	} else {
	    sharedMeshShape = mesh_from_internal(&validInternal.local, source);
	    if (sharedMeshShape)
		assign_shared_geometry_identity(sharedMeshShape,
		    dp->d_namep, typeLabel, revision, "surface");
	}

	if (!sharedVListShape && !sharedMeshShape) {
	    SbString msg;
	    msg.sprintf(
		"%s: direct leaf mesh conversion/tessellation failed for primitive type '%s'",
		fullPath.c_str(), typeLabel ? typeLabel : "");
	    source->realizationDiagnostic = msg;
	    return -1;
	}

	if (sharedVListShape) {
	    cache->storeMeshVListGeometry(cacheKey, sharedVListShape);
	    typeLabel = sharedVListShape->sourceType.getValue().getString();
	} else {
	    cache->storeMeshGeometry(cacheKey, sharedMeshShape);
	    typeLabel = sharedMeshShape->sourceType.getValue().getString();
	}
    }

    SoSeparator *leaf = realize_matrix_leaf_separator(pathMatrix);
    if (sharedVListShape) {
	SoBRLVListShape *vlistShape = new SoBRLVListShape;
	assign_realized_identity(vlistShape, NULL, fullPath.c_str(),
				 dp->d_namep, typeLabel, revision, source);
	vlistShape->setSharedGeometry(sharedVListShape);
	const char *geometryKind =
	    sharedVListShape->geometryKind.getValue().getString();
	vlistShape->geometryKind = (geometryKind && geometryKind[0]) ?
				   geometryKind : "line";
	if (geometryKind && BU_STR_EQUAL(geometryKind, "annotation"))
	    vlistShape->sourceType = "annotation";
	leaf->addChild(vlistShape);
    } else {
	SoBRLMeshShape *shape = mesh_instance_for_shared_geometry(
				    sharedMeshShape);
	assign_realized_identity(shape, NULL, fullPath.c_str(),
				 dp->d_namep, typeLabel, revision, source);
	shape->setSharedGeometry(sharedMeshShape);
	const char *geometryKind =
	    sharedMeshShape->geometryKind.getValue().getString();
	shape->geometryKind = (geometryKind && geometryKind[0]) ?
			      geometryKind : "surface";
	if (typeLabel && BU_STR_EQUAL(typeLabel, "bot") &&
	    source->lodBotThreshold.getValue() > 0)
	    publish_lod_mesh_if_available(shape, source, dbip, dp->d_namep);
	leaf->addChild(shape);
    }
    database_source_add_realized_child(source, leaf);
    return 1;
}

static int
realize_direct_leaf_mesh_compact(
    SoBRLDatabaseSource *source,
    BObolDatabaseSourceRealizationCache *cache,
    uint32_t revision,
    BObolCompactOccurrenceStream *stream)
{
    struct db_i *dbip = source ? source->getDatabase() : NULL;
    if (!source || !cache || !dbip || source_has_auxiliary_children(source))
	return 0;

    const std::string leafName =
	database_source_leaf_component(source->path.getValue());
    if (leafName.empty())
	return 0;

    struct directory *dp =
	db_lookup(dbip, leafName.c_str(), LOOKUP_QUIET);
    if (!dp || (dp->d_flags & RT_DIR_COMB))
	return 0;
    /* A streamed BREP needs its box-first, bounded producer path.  The direct
     * optimization would otherwise complete tessellation and PoP generation
     * before publishing any useful visual. */
    if (stream && source->lodBotThreshold.getValue() > 0 &&
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP)
	return 0;

    const std::string fullPath =
	database_source_full_path_string(source->path.getValue());
    SbMatrix pathMatrix;
    if (!database_source_path_matrix(dbip, source->path.getValue(),
	pathMatrix))
	return 0;
    SoBRLVListShape *sharedVListShape = NULL;
    SoBRLMeshShape *sharedMeshShape = NULL;
    SbBox3f cacheBounds;
    std::string cacheKey = realize_geometry_cache_key(dp);
    source_lod_cache_key_append(cacheKey, source, cacheBounds,
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP);
    BObolRealizationCacheMap<SoBRLVListShape *>::iterator foundVList =
	cache->sharedMeshVListGeometry.find(cacheKey);
    if (foundVList != cache->sharedMeshVListGeometry.end()) {
	sharedVListShape = foundVList->second;
    } else {
	BObolRealizationCacheMap<SoBRLMeshShape *>::iterator foundMesh =
	    cache->sharedMeshGeometry.find(cacheKey);
	if (foundMesh != cache->sharedMeshGeometry.end())
	    sharedMeshShape = foundMesh->second;
    }
    const BObolCachedPartGeometry *cachedWire =
	cache->findMeshVListCadGeometry(cacheKey);
    const BObolCachedPartGeometry *cachedMesh =
	cache->findMeshCadGeometry(cacheKey);
    if (sharedVListShape &&
	!source_cached_wire_matches_mesh_presentation(source,
	    sharedVListShape->sourceType.getValue().getString(),
	    sharedVListShape->geometryKind.getValue().getString()))
	sharedVListShape = NULL;
    if (cachedWire && !source_cached_wire_matches_mesh_presentation(source,
	    cachedWire->sourceType.c_str(), cachedWire->geometryKind.c_str()))
	cachedWire = NULL;
    if (!source_cached_mesh_matches_presentation(source, dp)) {
	sharedMeshShape = NULL;
	cachedMesh = NULL;
    }

    /*
     * Normal compact leaves do not need a Coin mesh or vlist carrier.  Mesh
     * LoD metadata and its optional cached display payload are value-backed
     * below, so the threshold does not require a temporary shape.
     */
    const bool hadCachedCad = cachedWire || cachedMesh;
    bool generatedViewDependentCsgGeometry = false;
    if (!sharedVListShape && !sharedMeshShape && !hadCachedCad) {
	owned_leaf_internal validInternal;
	if (rt_db_get_internal(&validInternal.local, dp, dbip, NULL) < 0 ||
	    !internal_payload_magic_valid(&validInternal.local)) {
	    SbString msg;
	    msg.sprintf("%s: direct compact leaf mesh internal fetch failed",
		fullPath.c_str());
	    source->realizationDiagnostic = msg;
	    if (validInternal.local.idb_ptr)
		rt_db_free_internal(&validInternal.local);
	    return -1;
	}
	validInternal.ownsLocal = true;

	const int internalType = validInternal.local.idb_type;
	const char *directTypeLabel = primitive_type_label(&validInternal.local);
	if (internalType == ID_MATERIAL) {
	    SoBRLMaterialObject *materialObject = material_object_from_internal(
		static_cast<struct rt_material_internal *>(
		    validInternal.local.idb_ptr));
	    if (!materialObject) {
		SbString msg;
		msg.sprintf("%s: material object realization failed",
		    fullPath.c_str());
		source->realizationDiagnostic = msg;
		return -1;
	    }
	    SoSeparator *leaf = realize_matrix_leaf_separator(pathMatrix);
	    assign_material_identity(materialObject, fullPath.c_str(),
		dp->d_namep, directTypeLabel, revision);
	    leaf->addChild(materialObject);
	    database_source_add_realized_child(source, leaf);
	    source->clearSourceBounds();
	    return 1;
	}

	const int drawMode = source_record_draw_mode(source);
	const bool wireGeometry = primitive_uses_wire_in_mesh_mode(internalType) ||
	    ((drawMode == BOBOL_LOD_DRAW_WIRE ||
	      (drawMode == BOBOL_LOD_DRAW_SHADED_BOTS &&
	       (!validInternal.local.idb_meth ||
		!validInternal.local.idb_meth->ft_indexed_face_set))) &&
	     internalType != ID_BOT);
	{
	    SbBox3f localBounds;
	    if (source_view_lod_active(source))
		(void)local_bounds_from_internal(&validInternal.local, localBounds);
	    Obol::PartGeometryBuilder generated;
	    int generatedGeometry = 0;
	    BObolSourceMeshRequest sourceMeshRequest;
	    bool lodBacked = false;
	    if (internalType == ID_PNTS) {
		generatedGeometry = cad_points_part_geometry_from_pnts(
		    static_cast<const struct rt_pnts_internal *>(
			validInternal.local.idb_ptr), generated);
	    } else if (wireGeometry) {
		generatedGeometry = cad_wire_part_geometry_from_lod_realization_internal(
		    &validInternal.local, source, localBounds, generated,
		    &generatedViewDependentCsgGeometry);
		if (!generatedGeometry)
		    generatedGeometry = cad_wire_part_geometry_from_plot_internal(
			&validInternal.local, source, generated);
	    } else {
		const struct rt_bot_internal *bot = internalType == ID_BOT ?
		    static_cast<const struct rt_bot_internal *>(
			validInternal.local.idb_ptr) : NULL;
		if (internalType == ID_BREP &&
		    source->lodBotThreshold.getValue() > 0 &&
		    validInternal.local.idb_meth &&
		    validInternal.local.idb_meth->ft_indexed_face_set) {
		    const struct bg_tess_tol ttol = source_tess_tol(source);
		    const struct bn_tol tol = BN_TOL_INIT_TOL;
		    std::shared_ptr<BObolStagedSourceMesh> staged =
			cad_staged_mesh_from_primitive_face_set(
			    dbip, dp, &validInternal.local, &ttol, &tol,
			    revision, sourceMeshRequest);
		    if (staged) {
			sourceMeshRequest.meshAssetPath = fullPath.c_str();
			sourceMeshRequest.meshAssetName = dp->d_namep;
			struct BObolMeshLodCacheStatus cacheStatus =
			    BOBOL_MESH_LOD_CACHE_STATUS_INIT;
			if (bobol_mesh_lod_cache_store_mesh_variant(
				dbip, dp->d_namep,
				staged->points, staged->pointCount,
				staged->normals, staged->faces,
				staged->faceCount, staged->contentKey,
				staged->shadedCullBackfaces,
				&cacheStatus) == BRLCAD_OK &&
			    cad_wire_part_geometry_from_aabb(
				sourceMeshRequest.bounds, generated)) {
			    sourceMeshRequest.meshAssetContentHash =
				cacheStatus.cache_key;
			    generatedGeometry = 1;
			    lodBacked = true;
			}
		    }
		}
		if (!generatedGeometry && bot &&
		    source->lodBotThreshold.getValue() > 0 &&
		    bot->num_faces >= source->lodBotThreshold.getValue() &&
		    cad_source_mesh_request_from_bot(sourceMeshRequest, bot)) {
		    sourceMeshRequest.meshAssetPath = fullPath.c_str();
		    sourceMeshRequest.meshAssetName = dp->d_namep;
		    if (cad_wire_part_geometry_from_aabb(sourceMeshRequest.bounds,
			    generated)) {
			generatedGeometry = 1;
			lodBacked = true;
		    }
		}
		if (!generatedGeometry) {
		    generatedGeometry = cad_mesh_part_geometry_from_internal(
			&validInternal.local, source, generated);
	    }
		if (generatedGeometry &&
		    drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE)
		    (void)cad_mesh_append_hidden_line_edges(generated);
	    }
		if (generatedGeometry) {
		    const char *geometryKind = internalType == ID_PNTS ? "point" :
			(wireGeometry ? "line" :
			 (lodBacked ? "aabb" : "surface"));
		if (primitive_is_annotation(internalType, directTypeLabel)) {
		    directTypeLabel = "annotation";
		    geometryKind = "annotation";
		}
		if (wireGeometry || internalType == ID_PNTS)
		    cache->storeMeshVListCadGeometry(cacheKey,
			std::move(generated), directTypeLabel, geometryKind,
			localBounds.isEmpty() ? NULL : &localBounds, false, NULL,
			generatedViewDependentCsgGeometry);
		else
		    cache->storeMeshCadGeometry(cacheKey, std::move(generated),
			directTypeLabel, geometryKind, NULL, lodBacked,
			lodBacked ? &sourceMeshRequest : NULL);
		cachedWire = cache->findMeshVListCadGeometry(cacheKey);
		cachedMesh = cache->findMeshCadGeometry(cacheKey);
	    }
	}
    }
    bobol_performance_counter_add(
	(sharedVListShape || sharedMeshShape || hadCachedCad) ?
	BOBOL_PERF_MESH_CACHE_HITS :
	BOBOL_PERF_MESH_CACHE_MISSES, 1);

    const char *typeLabel = sharedVListShape ?
			    sharedVListShape->sourceType.getValue().getString() :
			    (sharedMeshShape ?
			     sharedMeshShape->sourceType.getValue().getString() :
			     (cachedWire && !cachedWire->sourceType.empty() ?
			      cachedWire->sourceType.c_str() :
			      (cachedMesh && !cachedMesh->sourceType.empty() ?
			       cachedMesh->sourceType.c_str() : NULL)));
    if (!sharedVListShape && !sharedMeshShape && !cachedWire && !cachedMesh) {
	owned_leaf_internal validInternal;
	if (rt_db_get_internal(&validInternal.local, dp, dbip, NULL) < 0 ||
	    !internal_payload_magic_valid(&validInternal.local)) {
	    SbString msg;
	    msg.sprintf("%s: direct compact leaf mesh internal fetch failed",
			fullPath.c_str());
	    source->realizationDiagnostic = msg;
	    if (validInternal.local.idb_ptr)
		rt_db_free_internal(&validInternal.local);
	    return -1;
	}
	validInternal.ownsLocal = true;

	typeLabel = primitive_type_label(&validInternal.local);
	const int internalType = validInternal.local.idb_type;
	SbBox3f localBounds;
	if (source_view_lod_active(source))
	    (void)local_bounds_from_internal(&validInternal.local, localBounds);
	if (internalType == ID_MATERIAL) {
	    SoBRLMaterialObject *materialObject = material_object_from_internal(
		static_cast<struct rt_material_internal *>(
		    validInternal.local.idb_ptr));
	    if (!materialObject) {
		SbString msg;
		msg.sprintf("%s: material object realization failed",
		    fullPath.c_str());
		source->realizationDiagnostic = msg;
		return -1;
	    }
	    SoSeparator *leaf = new SoSeparator;
	    assign_material_identity(materialObject, fullPath.c_str(),
		dp->d_namep, typeLabel, revision);
	    leaf->addChild(materialObject);
	    database_source_add_realized_child(source, leaf);
	    source->clearSourceBounds();
	    return 1;
	}

	if (internalType == ID_PNTS) {
	    sharedVListShape = vlist_from_pnts(
		static_cast<const struct rt_pnts_internal *>(
		    validInternal.local.idb_ptr));
	    if (sharedVListShape)
		assign_shared_geometry_identity(sharedVListShape,
		    dp->d_namep, typeLabel, revision, "point");
	} else if (primitive_uses_wire_in_mesh_mode(internalType)) {
	    sharedVListShape = vlist_from_lod_realization_internal(
		&validInternal.local, source, localBounds);
	    if (!sharedVListShape)
		sharedVListShape = vlist_from_plot_internal(&validInternal.local,
		    source);
	    if (sharedVListShape) {
		assign_shared_geometry_identity(sharedVListShape,
		    dp->d_namep, typeLabel, revision, "line");
		if (primitive_is_annotation(internalType, typeLabel)) {
		    sharedVListShape->sourceType = "annotation";
		    sharedVListShape->geometryKind = "annotation";
		}
	    }
	} else if ((source_record_draw_mode(source) == BOBOL_LOD_DRAW_WIRE ||
		    (source_record_draw_mode(source) ==
			BOBOL_LOD_DRAW_SHADED_BOTS &&
		     (!validInternal.local.idb_meth ||
		      !validInternal.local.idb_meth->ft_indexed_face_set))) &&
			   internalType != ID_BOT) {
	    sharedVListShape = vlist_from_lod_realization_internal(
		&validInternal.local, source, localBounds);
	    if (!sharedVListShape)
		sharedVListShape = vlist_from_plot_internal(&validInternal.local,
		    source);
	    if (sharedVListShape)
		assign_shared_geometry_identity(sharedVListShape,
		    dp->d_namep, typeLabel, revision, "line");
	} else {
	    sharedMeshShape = mesh_from_internal(&validInternal.local, source);
	    if (sharedMeshShape)
		assign_shared_geometry_identity(sharedMeshShape,
		    dp->d_namep, typeLabel, revision, "surface");
	}

	if (!sharedVListShape && !sharedMeshShape) {
	    SbString msg;
	    msg.sprintf(
		"%s: direct compact leaf mesh conversion/tessellation failed for primitive type '%s'",
		fullPath.c_str(), typeLabel ? typeLabel : "");
	    source->realizationDiagnostic = msg;
	    return -1;
	}

	if (sharedVListShape) {
	    typeLabel = sharedVListShape->sourceType.getValue().getString();
	} else {
	    typeLabel = sharedMeshShape->sourceType.getValue().getString();
	}
    }

    BObolCompactOccurrence occurrence;
    if (sharedVListShape || cachedWire) {
	const char *geometryKind = sharedVListShape ?
	    sharedVListShape->geometryKind.getValue().getString() :
	    (cachedWire && !cachedWire->geometryKind.empty() ?
	     cachedWire->geometryKind.c_str() : "line");
	occurrence.geometry = cachedWire ? cachedWire->geometry :
	    std::shared_ptr<const Obol::PartGeometry>();
	occurrence.viewDependentCsgGeometry = cachedWire ?
	    (cachedWire->viewDependentCsgGeometry ? TRUE : FALSE) :
	    (generatedViewDependentCsgGeometry ? TRUE : FALSE);
	if (!occurrence.geometry && sharedVListShape) {
	    Obol::PartGeometryBuilder generated;
	    if (cad_vlist_part_geometry(sharedVListShape, generated))
		occurrence.geometry = cache->storeMeshVListCadGeometry(cacheKey,
		    std::move(generated), typeLabel, geometryKind);
	}
	occurrence.summary = compact_occurrence_summary(source,
	    fullPath.c_str(), dp->d_namep,
	    geometryKind && BU_STR_EQUAL(geometryKind, "annotation") ?
	    "annotation" : typeLabel,
	    geometryKind, revision, BObolRealizedShapeSummary::SHAPE_VLIST);
    } else {
	const char *geometryKind = sharedMeshShape ?
	    sharedMeshShape->geometryKind.getValue().getString() :
	    (cachedMesh && !cachedMesh->geometryKind.empty() ?
	     cachedMesh->geometryKind.c_str() : "surface");
	if (sharedMeshShape) {
	    sharedMeshShape->drawMode = source_record_draw_mode(source);
	    sharedMeshShape->hiddenLine =
		source_record_draw_mode(source) ==
		BOBOL_LOD_DRAW_HIDDEN_LINE ? TRUE : FALSE;
	    if (typeLabel && BU_STR_EQUAL(typeLabel, "bot") &&
		source->lodBotThreshold.getValue() > 0)
		publish_lod_mesh_if_available(sharedMeshShape, source, dbip,
		    dp->d_namep);
	}
	occurrence.geometry = cachedMesh ? cachedMesh->geometry :
	    std::shared_ptr<const Obol::PartGeometry>();
	if (!occurrence.geometry && sharedMeshShape) {
	    Obol::PartGeometryBuilder generated;
	    if (cad_mesh_part_geometry(sharedMeshShape, generated))
		occurrence.geometry = cache->storeMeshCadGeometry(cacheKey,
		    std::move(generated), typeLabel, geometryKind, NULL,
		    sharedMeshShape->isLodBackedMesh());
	}
	occurrence.lodBacked = sharedMeshShape ?
	    sharedMeshShape->isLodBackedMesh() :
	    (cachedMesh && cachedMesh->lodBacked ? TRUE : FALSE);
	if (sharedMeshShape)
	    occurrence.sourceMeshRequestValid =
		sharedMeshShape->makeSourceMeshRequest(
		    occurrence.sourceMeshRequest);
	else if (cachedMesh && cachedMesh->sourceMeshRequestValid) {
	    occurrence.sourceMeshRequestValid = TRUE;
	    occurrence.sourceMeshRequest = cachedMesh->sourceMeshRequest;
	}
	occurrence.summary = compact_occurrence_summary(source,
	    fullPath.c_str(), dp->d_namep, typeLabel, geometryKind, revision,
	    BObolRealizedShapeSummary::SHAPE_MESH);
	if (sharedMeshShape) {
	    BObolRealizedShapeSummary meshSummary;
	    realized_mesh_shape_summary(sharedMeshShape, meshSummary);
	    occurrence.summary.lodAvailable = meshSummary.lodAvailable;
	    occurrence.summary.lodActiveCut = meshSummary.lodActiveCut;
	    occurrence.summary.lodFaceCount = meshSummary.lodFaceCount;
	    occurrence.summary.lodPointCount = meshSummary.lodPointCount;
	    occurrence.summary.lodOriginalPointCount =
		meshSummary.lodOriginalPointCount;
	    occurrence.summary.lodNormalCount = meshSummary.lodNormalCount;
	    occurrence.summary.lodHasSnappedPoints =
		meshSummary.lodHasSnappedPoints;
	    occurrence.summary.lodHasNormals = meshSummary.lodHasNormals;
	    occurrence.summary.lodBoundsMin = meshSummary.lodBoundsMin;
	    occurrence.summary.lodBoundsMax = meshSummary.lodBoundsMax;
	}
	if (occurrence.sourceMeshRequestValid)
	    compact_source_mesh_request_sync(occurrence.sourceMeshRequest,
		occurrence.summary);
	if (occurrence.sourceMeshRequestValid)
	    compact_summary_lod_from_source_mesh_request(occurrence.summary,
		occurrence.sourceMeshRequest);
	if (occurrence.sourceMeshRequestValid)
	    cache_mesh_cad_source_request(cache, cacheKey,
		occurrence.sourceMeshRequest);
    }

    occurrence.occurrenceIndex = source->occurrenceIndex.getValue();
    occurrence.booleanOperation = source->booleanOperation.getValue();
    occurrence.localTransform = pathMatrix;
    const int compacted = source->setCompactOccurrence(occurrence);
    if (compacted > 0)
	publish_direct_compact_occurrence(source, stream, dp, occurrence);
    return compacted > 0 ? 1 : 0;
}

static bool
compact_stream_pca_bucket(std::array<int64_t, 3> &bucket,
	const struct bg_pca_frame &frame)
{
    for (size_t axis = 0; axis < bucket.size(); axis++) {
	const double value = static_cast<double>(frame.singular_values[axis]);
	if (value <= SMALL_FASTF) {
	    bucket[axis] = std::numeric_limits<int64_t>::min();
	    continue;
	}
	const double scaled = log2(value) * 4096.0;
	if (!std::isfinite(scaled) ||
	    scaled > static_cast<double>(std::numeric_limits<int64_t>::max()) ||
	    scaled < static_cast<double>(std::numeric_limits<int64_t>::min()))
	    return false;
	bucket[axis] = static_cast<int64_t>(llround(scaled));
    }
    return true;
}

static bool
compact_stream_lod_pca_signature(
	compact_stream_lod_reuse_entry &entry,
	const struct rt_bot_internal *bot)
{
    if (entry.signatureValid)
	return true;
    if (!bot || !bot->vertices || !bot->faces ||
	!bot->num_vertices || !bot->num_faces ||
	bg_trimesh_pca_get_signature(&entry.signature, bot->faces,
	    bot->num_faces, reinterpret_cast<const point_t *>(bot->vertices),
	    bot->num_vertices, VUNITIZE_TOL, 1.0e-6) != BRLCAD_OK ||
	!compact_stream_pca_bucket(entry.bucket, entry.signature.frame))
	return false;

    const double characteristicExtent =
	static_cast<double>(entry.signature.frame.singular_values[0]) /
	sqrt(static_cast<double>(bot->num_vertices));
    if (std::isfinite(characteristicExtent) && characteristicExtent > 0.0)
	entry.comparisonTolerance = static_cast<fastf_t>(std::max(
	    static_cast<double>(VUNITIZE_TOL),
	    characteristicExtent * 1.0e-9));
    entry.signatureValid = true;
    return true;
}

/* Exact count matching alone is a poor PCA broad phase for generated vehicle
 * meshes: thousands of unrelated BoTs commonly share one topology size.  Hash
 * a bounded, deterministic sample of edge lengths before scanning the full
 * mesh.  Squared edge lengths are invariant under the rigid transforms PCA
 * reuse accepts, and quantized logarithms tolerate the roundoff introduced by
 * xpush.  A match is only permission to run the existing exact PCA/equality
 * check, so collisions cannot cause incorrect geometry reuse.  Vertex/face
 * reordering may produce a safe false negative and simply forgo optimization. */
static bool
compact_stream_lod_sample_fingerprint(uint64_t &fingerprint,
	const struct rt_bot_internal *bot)
{
    fingerprint = 0;
    if (!bot || !bot->vertices || !bot->faces ||
	!bot->num_vertices || !bot->num_faces)
	return false;

    const size_t sampleCount = std::min<size_t>(32, bot->num_faces);
    uint64_t hash = 1469598103934665603ULL;
    const auto mix = [&hash](uint64_t word) {
	hash ^= word;
	hash *= 1099511628211ULL;
    };
    mix(static_cast<uint64_t>(bot->num_vertices));
    mix(static_cast<uint64_t>(bot->num_faces));
    for (size_t sample = 0; sample < sampleCount; sample++) {
	const size_t face = sampleCount > 1 ?
	    sample * (bot->num_faces - 1) / (sampleCount - 1) : 0;
	const int *indices = &bot->faces[face * 3];
	for (size_t edge = 0; edge < 3; edge++) {
	    const int first = indices[edge];
	    const int second = indices[(edge + 1) % 3];
	    if (first < 0 || second < 0 ||
		static_cast<size_t>(first) >= bot->num_vertices ||
		static_cast<size_t>(second) >= bot->num_vertices)
		return false;
	    const fastf_t *a = &bot->vertices[static_cast<size_t>(first) * 3];
	    const fastf_t *b = &bot->vertices[static_cast<size_t>(second) * 3];
	    const double dx = static_cast<double>(a[X] - b[X]);
	    const double dy = static_cast<double>(a[Y] - b[Y]);
	    const double dz = static_cast<double>(a[Z] - b[Z]);
	    const double lengthSquared = dx * dx + dy * dy + dz * dz;
	    if (!std::isfinite(lengthSquared) || lengthSquared < 0.0)
		return false;
	    int64_t quantized = std::numeric_limits<int64_t>::min();
	    if (lengthSquared > SMALL_FASTF) {
		const double scaled = log2(lengthSquared) * 262144.0;
		if (!std::isfinite(scaled) ||
		    scaled > static_cast<double>(
			std::numeric_limits<int64_t>::max()) ||
		    scaled < static_cast<double>(
			std::numeric_limits<int64_t>::min()))
		    return false;
		quantized = static_cast<int64_t>(llround(scaled));
	    }
	    mix(static_cast<uint64_t>(quantized));
	}
    }
    fingerprint = hash;
    return true;
}

static void
compact_stream_lod_asset_record_store(struct db_i *dbip,
	struct directory *objectDp,
	const BObolSourceMeshRequest &objectRequest,
	const BObolSourceMeshRequest &assetRequest,
	const SbMatrix &assetToObject)
{
    if (!dbip || !objectDp || !objectDp->d_namep ||
	assetRequest.meshAssetName.getLength() == 0 ||
	objectRequest.bounds.isEmpty() || assetRequest.meshAssetBounds.isEmpty())
	return;

    BObolDrawLodAssetRecord record;
    bobol_draw_lod_asset_record_init(&record);
    bu_strlcpy(record.assetName,
	assetRequest.meshAssetName.getString(), sizeof(record.assetName));
    record.faceCount = assetRequest.faceCount;
    record.pointCount = assetRequest.pointCount;
    const SbVec3f objectMin = objectRequest.bounds.getMin();
    const SbVec3f objectMax = objectRequest.bounds.getMax();
    const SbVec3f assetMin = assetRequest.meshAssetBounds.getMin();
    const SbVec3f assetMax = assetRequest.meshAssetBounds.getMax();
    VSET(record.boundsMin, objectMin[0], objectMin[1], objectMin[2]);
    VSET(record.boundsMax, objectMax[0], objectMax[1], objectMax[2]);
    VSET(record.assetBoundsMin, assetMin[0], assetMin[1], assetMin[2]);
    VSET(record.assetBoundsMax, assetMax[0], assetMax[1], assetMax[2]);
    for (size_t row = 0; row < 4; row++)
	for (size_t column = 0; column < 4; column++)
	    record.assetToObject[column * 4 + row] =
		static_cast<fastf_t>(assetToObject[row][column]);
    (void)bobol_draw_lod_asset_cache_store(dbip, objectDp->d_namep,
	&record);
}

static const BObolCachedPartGeometry *
compact_stream_lod_cached_asset(realize_walk_data *data,
	struct db_i *dbip, struct directory *dp,
	const std::string &objectCacheKey)
{
    if (!data || !data->cache || !data->source || !dbip || !dp ||
	!dp->d_namep || data->source->lodBotThreshold.getValue() == 0)
	return NULL;

    BObolDrawLodAssetRecord record;
    if (bobol_draw_lod_asset_cache_get(dbip, dp->d_namep, &record) !=
	BRLCAD_OK ||
	record.faceCount < data->source->lodBotThreshold.getValue())
	return NULL;

    struct directory *assetDp =
	db_lookup(dbip, record.assetName, LOOKUP_QUIET);
    if (assetDp == RT_DIR_NULL)
	return NULL;
    BObolSourceMeshRequest assetRequest;
    assetRequest.faceCount = record.faceCount;
    assetRequest.pointCount = record.pointCount;
    assetRequest.bounds = SbBox3f(
	SbVec3f(record.assetBoundsMin[X], record.assetBoundsMin[Y],
	    record.assetBoundsMin[Z]),
	SbVec3f(record.assetBoundsMax[X], record.assetBoundsMax[Y],
	    record.assetBoundsMax[Z]));
    assetRequest.meshAssetBounds = assetRequest.bounds;
    assetRequest.meshAssetName = record.assetName;
    assetRequest.meshAssetPath = record.assetName;

    std::string assetCacheKey = realize_geometry_cache_key(assetDp);
    SbBox3f unusedBounds;
    source_lod_cache_key_append(assetCacheKey, data->source, unusedBounds);
    const BObolCachedPartGeometry *assetGeometry =
	data->cache->findMeshCadGeometry(assetCacheKey);
    if (!assetGeometry) {
	Obol::PartGeometryBuilder generated;
	if (!cad_wire_part_geometry_from_aabb(assetRequest.bounds, generated))
	    return NULL;
	data->cache->storeMeshCadGeometry(assetCacheKey, std::move(generated),
	    "bot", "aabb", &assetRequest.bounds, true, &assetRequest);
	assetGeometry = data->cache->findMeshCadGeometry(assetCacheKey);
    }
    if (!assetGeometry || !assetGeometry->geometry ||
	!assetGeometry->sourceMeshRequestValid)
	return NULL;

    const SbMatrix assetToObject = mat_to_sbmatrix(record.assetToObject);
    if (objectCacheKey != assetCacheKey) {
	const SbBox3f objectBounds(
	    SbVec3f(record.boundsMin[X], record.boundsMin[Y],
		record.boundsMin[Z]),
	    SbVec3f(record.boundsMax[X], record.boundsMax[Y],
		record.boundsMax[Z]));
	BObolSourceMeshRequest objectRequest =
	    assetGeometry->sourceMeshRequest;
	objectRequest.meshAssetTransform = assetToObject;
	data->cache->storeMeshCadGeometryReference(objectCacheKey,
	    assetGeometry->geometry, assetToObject, "bot", "aabb",
	    &objectBounds, true, &objectRequest);
    }
    data->stream_lod_asset_hits++;
    if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
	bu_log("[obol-timing] stream LoD asset cache: %s -> %s\n",
	    dp->d_namep, record.assetName);
    return data->cache->findMeshCadGeometry(objectCacheKey);
}

static bool
compact_stream_lod_transformed_reuse(
	realize_walk_data *data, struct db_i *dbip, struct directory *dp,
	const char *path, const std::string &cacheKey,
	const struct rt_bot_internal *bot,
	const BObolSourceMeshRequest &sourceMeshRequest)
{
    if (!data || !data->cache || !dbip || !dp || !bot || !bot->vertices ||
	!bot->faces || !bot->num_vertices || !bot->num_faces ||
	(bot->bot_flags & RT_BOT_HAS_SURFACE_NORMALS))
	return false;

    compact_stream_lod_reuse_entry candidate;
    candidate.dp = dp;
    candidate.cacheKey = cacheKey;
    candidate.sourceMeshRequest = sourceMeshRequest;
    candidate.sourceMeshRequest.meshAssetPath = path ? path : "";
    candidate.sourceMeshRequest.meshAssetName =
	dp->d_namep ? dp->d_namep : "";
    candidate.vertexCount = bot->num_vertices;
    candidate.faceCount = bot->num_faces;
    candidate.mode = bot->mode;
    candidate.orientation = bot->orientation;
    candidate.sampleFingerprintValid =
	compact_stream_lod_sample_fingerprint(candidate.sampleFingerprint, bot);

    /* Face/vertex counts, BoT mode, and orientation are exact, essentially
     * free broad-phase invariants.  The bounded edge fingerprint then
     * distinguishes the common case where thousands of unrelated meshes share
     * those counts.  Do not make those shapes pay for a full PCA traversal:
     * keep the first as an unevaluated representative and calculate signatures
     * only when a later occurrence could actually be a transformed copy. */
    bool haveBroadCandidate = false;
    for (const compact_stream_lod_reuse_entry &representative :
	data->stream_lod_reuse) {
	if (representative.vertexCount == candidate.vertexCount &&
	    representative.faceCount == candidate.faceCount &&
	    representative.mode == candidate.mode &&
	    representative.orientation == candidate.orientation &&
	    (!candidate.sampleFingerprintValid ||
	     !representative.sampleFingerprintValid ||
	     representative.sampleFingerprint == candidate.sampleFingerprint)) {
	    haveBroadCandidate = true;
	    break;
	}
    }
    if (!haveBroadCandidate) {
	compact_stream_lod_asset_record_store(dbip, dp, sourceMeshRequest,
	    candidate.sourceMeshRequest, SbMatrix::identity());
	data->stream_lod_pca_deferred++;
	if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
	    bu_log("[obol-timing] stream PCA: %s registered deferred "
		   "(no broad candidate)\n",
		   dp->d_namep ? dp->d_namep : "?");
	data->stream_lod_reuse.push_back(std::move(candidate));
	return false;
    }

    if (!compact_stream_lod_pca_signature(candidate, bot)) {
	compact_stream_lod_asset_record_store(dbip, dp, sourceMeshRequest,
	    candidate.sourceMeshRequest, SbMatrix::identity());
	data->stream_lod_reuse.push_back(std::move(candidate));
	return false;
    }
    data->stream_lod_pca_evaluated++;

    size_t plausibleCandidates = 0;
    for (compact_stream_lod_reuse_entry &representative :
	data->stream_lod_reuse) {
	if (representative.vertexCount != candidate.vertexCount ||
	    representative.faceCount != candidate.faceCount ||
	    representative.mode != candidate.mode ||
	    representative.orientation != candidate.orientation ||
	    (candidate.sampleFingerprintValid &&
	     representative.sampleFingerprintValid &&
	     representative.sampleFingerprint != candidate.sampleFingerprint))
	    continue;

	if (!data->stream_lod_cached_representative_valid ||
	    data->stream_lod_cached_representative_dp != representative.dp) {
	    if (data->stream_lod_cached_representative_valid) {
		rt_db_free_internal(
		    &data->stream_lod_cached_representative);
		RT_DB_INTERNAL_INIT(
		    &data->stream_lod_cached_representative);
		data->stream_lod_cached_representative_valid = false;
		data->stream_lod_cached_representative_dp = NULL;
	    }
	    if (rt_db_get_internal(
		    &data->stream_lod_cached_representative,
		    representative.dp, dbip, NULL) < 0 ||
		data->stream_lod_cached_representative.idb_type != ID_BOT ||
		!data->stream_lod_cached_representative.idb_ptr) {
		if (data->stream_lod_cached_representative.idb_ptr)
		    rt_db_free_internal(
			&data->stream_lod_cached_representative);
		RT_DB_INTERNAL_INIT(
		    &data->stream_lod_cached_representative);
		continue;
	    }
	    data->stream_lod_cached_representative_dp = representative.dp;
	    data->stream_lod_cached_representative_valid = true;
	    data->stream_lod_representative_imports++;
	}
	const struct rt_bot_internal *representativeBot =
	    static_cast<const struct rt_bot_internal *>(
		data->stream_lod_cached_representative.idb_ptr);
	if (!representativeBot) {
	    continue;
	}
	if (!compact_stream_lod_pca_signature(representative,
		representativeBot)) {
	    continue;
	}
	bool nearby = true;
	for (size_t axis = 0; axis < candidate.bucket.size(); axis++) {
	    if (representative.bucket[axis] ==
		    std::numeric_limits<int64_t>::min() ||
		candidate.bucket[axis] ==
		    std::numeric_limits<int64_t>::min()) {
		if (representative.bucket[axis] != candidate.bucket[axis])
		    nearby = false;
		continue;
	    }
	    if (llabs(representative.bucket[axis] -
		    candidate.bucket[axis]) > 1)
		nearby = false;
	}
	if (!nearby) {
	    continue;
	}
	plausibleCandidates++;
	const fastf_t tolerance = std::max(
	    representative.comparisonTolerance,
	    candidate.comparisonTolerance);
	const bool equal = bg_trimesh_pca_equal(&representative.signature,
	    representativeBot->faces, representativeBot->num_faces,
	    reinterpret_cast<const point_t *>(representativeBot->vertices),
	    representativeBot->num_vertices, &candidate.signature, bot->faces,
	    bot->num_faces, reinterpret_cast<const point_t *>(bot->vertices),
	    bot->num_vertices, tolerance) == 0;
	if (!equal)
	    continue;

	mat_t representativeToCandidate;
	if (bg_pca_frame_relative_matrix(representativeToCandidate,
		&representative.signature.frame,
		&candidate.signature.frame) != BRLCAD_OK)
	    continue;
	const BObolCachedPartGeometry *geometry =
	    data->cache->findMeshCadGeometry(representative.cacheKey);
	if (!geometry || !geometry->geometry ||
	    !geometry->sourceMeshRequestValid)
	    continue;
	const SbMatrix transform =
	    mat_to_sbmatrix(representativeToCandidate);
	compact_stream_lod_asset_record_store(dbip, dp, sourceMeshRequest,
	    geometry->sourceMeshRequest, transform);
	const SbBox3f bounds = database_source_transform_bounds(
	    geometry->bounds, transform);
	BObolSourceMeshRequest reusedRequest = geometry->sourceMeshRequest;
	reusedRequest.meshAssetTransform = transform;
	data->cache->storeMeshCadGeometryReference(cacheKey,
	    geometry->geometry, transform,
	    geometry->sourceType.empty() ? "bot" :
		geometry->sourceType.c_str(),
	    geometry->geometryKind.empty() ? "surface" :
		geometry->geometryKind.c_str(),
	    &bounds, true, &reusedRequest);
	data->stream_lod_pca_reused++;
	if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
	    bu_log("[obol-timing] stream PCA: %s reuses %s "
		   "(tol %.6g, %zu candidates)\n",
		   dp->d_namep ? dp->d_namep : "?",
		   representative.dp && representative.dp->d_namep ?
		       representative.dp->d_namep : "?",
		   tolerance, plausibleCandidates);
	return true;
    }

    compact_stream_lod_asset_record_store(dbip, dp, sourceMeshRequest,
	candidate.sourceMeshRequest, SbMatrix::identity());
    if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
	bu_log("[obol-timing] stream PCA: %s registered "
	       "(tol %.6g, %zu candidates)\n",
	       dp->d_namep ? dp->d_namep : "?",
	       candidate.comparisonTolerance, plausibleCandidates);
    data->stream_lod_reuse.push_back(std::move(candidate));
    return false;
}

static union tree *
    realize_mesh_leaf(struct db_tree_state *tsp,
		      const struct db_full_path *pathp,
		      struct directory *dp,
		      void *client_data)
{
    struct realize_walk_data *data = static_cast<struct realize_walk_data *>(client_data);
    if (!data || !data->source || !pathp || !tsp || !tsp->ts_dbip ||
	!dp)
	return TREE_NULL;

    /* A streaming realization that has been cancelled stops tessellating the
     * remaining leaves promptly; the partial result is discarded by the pump. */
    if (data->stream_sink && data->stream_sink->isCancelled())
	return make_nop_tree();

    std::string walkOccurrenceIdentity;
    uint32_t duplicateOrdinal = 0;
    if (data->compact_index) {
	const std::string identity = realize_walk_instance_identity(tsp, pathp);
	if (identity.empty()) {
	    data->compact_unsupported = 1;
	    return TREE_NULL;
	}
	if (!data->compact_seen_instances.insert(identity).second)
	    return make_nop_tree();
	walkOccurrenceIdentity = realize_walk_occurrence_identity(pathp);
	duplicateOrdinal = data->compact_occurrence_counts[
	    walkOccurrenceIdentity]++;
    }

    /*
     * Mesh prefill covers unthresholded ordinary shaded leaves.  Publish
     * point, wire, and thresholded mesh cases directly instead of making a
     * transient Coin carrier only to convert it into cached PartGeometry.
     */
    if (data->compact_index) {
	SbBox3f cacheBounds;
	std::string cacheKey = realize_geometry_cache_key(dp);
	source_lod_cache_key_append(cacheKey, data->source, cacheBounds,
	    dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP);
	const BObolCachedPartGeometry *cachedWire =
	    data->cache->findMeshVListCadGeometry(cacheKey);
	const BObolCachedPartGeometry *cachedMesh =
	    data->cache->findMeshCadGeometry(cacheKey);
	bool generatedViewDependentCsgGeometry = cachedWire ?
	    cachedWire->viewDependentCsgGeometry : false;
	if (cachedWire && !source_cached_wire_matches_mesh_presentation(
		data->source, cachedWire->sourceType.c_str(),
		cachedWire->geometryKind.c_str()))
	    cachedWire = NULL;
	if (!source_cached_mesh_matches_presentation(data->source, dp))
	    cachedMesh = NULL;
	const bool hadCachedWire = cachedWire != NULL;
	if (!cachedWire && !cachedMesh) {
	    if (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT)
		cachedMesh = compact_stream_lod_cached_asset(data,
		    tsp->ts_dbip, dp, cacheKey);
	}
	if (!cachedWire && !cachedMesh) {
	    owned_leaf_internal validInternal;
	    struct rt_db_internal *localIntern =
		import_walk_leaf_internal(tsp, dp, &validInternal);
	    if (!localIntern) {
		data->failed_shapes++;
		set_leaf_import_diagnostic(data, pathp,
			validInternal.ownsLocal ? &validInternal.local : NULL);
		return TREE_NULL;
	    }
	    if (database_source_terminal_empty_bot(localIntern))
		return make_nop_tree();

	    const int internalType = localIntern->idb_type;
	    const int drawMode = source_record_draw_mode(data->source);
	    const bool wireGeometry =
		primitive_uses_wire_in_mesh_mode(internalType) ||
		((drawMode == BOBOL_LOD_DRAW_WIRE ||
		  (drawMode == BOBOL_LOD_DRAW_SHADED_BOTS &&
		   (!localIntern->idb_meth ||
		    !localIntern->idb_meth->ft_indexed_face_set))) &&
		 internalType != ID_BOT);
	    if (internalType == ID_PNTS || wireGeometry) {
		SbBox3f localBounds;
		if (source_view_lod_active(data->source))
		    (void)local_bounds_from_internal(localIntern, localBounds);
		Obol::PartGeometryBuilder generated;
		int generatedGeometry = 0;
		if (internalType == ID_PNTS) {
		    generatedGeometry = cad_points_part_geometry_from_pnts(
			static_cast<const struct rt_pnts_internal *>(
			    localIntern->idb_ptr), generated);
		} else {
		    generatedGeometry =
			cad_wire_part_geometry_from_lod_realization_internal(
			    localIntern, data->source, localBounds, generated,
			    &generatedViewDependentCsgGeometry);
		    if (!generatedGeometry)
			generatedGeometry =
			    cad_wire_part_geometry_from_plot_internal(localIntern,
				data->source, generated);
		}
		if (generatedGeometry) {
		    const char *typeLabel = primitive_type_label(localIntern);
		    const char *geometryKind = internalType == ID_PNTS ? "point" :
			"line";
		    if (primitive_is_annotation(internalType, typeLabel)) {
			typeLabel = "annotation";
			geometryKind = "annotation";
		    }
		    data->cache->storeMeshVListCadGeometry(cacheKey,
			std::move(generated), typeLabel, geometryKind,
			localBounds.isEmpty() ? NULL : &localBounds, false, NULL,
			generatedViewDependentCsgGeometry);
		    cachedWire = data->cache->findMeshVListCadGeometry(cacheKey);
		}
	    } else if (data->source->lodBotThreshold.getValue() > 0) {
		Obol::PartGeometryBuilder generated;
		BObolSourceMeshRequest sourceMeshRequest;
		bool lodBacked = false;
		bool transformedReuse = false;
		const struct rt_bot_internal *bot = internalType == ID_BOT ?
		    static_cast<const struct rt_bot_internal *>(
			localIntern->idb_ptr) : NULL;
		int generatedGeometry = 0;
		if (bot && bot->num_faces >=
		    data->source->lodBotThreshold.getValue() &&
		    cad_source_mesh_request_from_bot(sourceMeshRequest, bot) &&
		    cad_wire_part_geometry_from_aabb(sourceMeshRequest.bounds,
			generated)) {
		    generatedGeometry = 1;
		    lodBacked = true;
		    char *assetPath = db_path_to_string(pathp);
		    sourceMeshRequest.meshAssetPath =
			assetPath ? assetPath : dp->d_namep;
		    sourceMeshRequest.meshAssetName = dp->d_namep;
		    /*
		     * The AABB is already useful drawing data.  Publish it before
		     * PCA matching/cache registration, which can be expensive for
		     * a first-seen very large BoT.  A later occurrence with the
		     * same path and AABB geometry upgrades this live entry in
		     * place by adding its source-mesh request; the LoD provider can
		     * then replace the box with progressively richer PoP data.
		     */
		    if (data->stream_sink &&
			!data->stream_sink->isCancelled()) {
			BObolCompactOccurrence provisional;
			provisional.geometry = bobol_cad_build_geometry(
			    generated, "provisional large-mesh bounds");
			if (provisional.geometry) {
			    provisional.localTransform =
				mat_to_sbmatrix(tsp->ts_mat);
			    std::string provisionalPath =
				assetPath ? assetPath :
				(dp->d_namep ? dp->d_namep : "");
			    if (duplicateOrdinal > 0) {
				char suffix[32] = {0};
				snprintf(suffix, sizeof(suffix), "@%u",
				    duplicateOrdinal);
				provisionalPath += suffix;
			    }
			    provisional.summary =
				compact_occurrence_tree_summary(
				    data->source, tsp, pathp,
				    provisionalPath.c_str(),
				    dp->d_namep,
				    primitive_type_label(localIntern),
				    "aabb", data->revision,
				    BObolRealizedShapeSummary::SHAPE_MESH,
				    static_cast<BObolMaterialColorSweep *>(
					data->material_sweep));
			    provisional.occurrenceIndex =
				pathp->fp_cinst && pathp->fp_len ?
				static_cast<uint32_t>(
				    DB_FULL_PATH_GET_COMB_INST(pathp,
					pathp->fp_len - 1)) : 0;
			    provisional.booleanOperation =
				(tsp->ts_sofar & TS_SOFAR_MINUS) ?
				SoBRLDatabaseSource::BOOLEAN_SUBTRACT :
				((tsp->ts_sofar & TS_SOFAR_INTER) ?
				 SoBRLDatabaseSource::BOOLEAN_INTERSECT :
				 SoBRLDatabaseSource::BOOLEAN_UNION);
			    realize_walk_stream_push(data, provisional);
			}
		    }
		    transformedReuse =
			compact_stream_lod_transformed_reuse(data,
			    tsp->ts_dbip, dp,
			    sourceMeshRequest.meshAssetPath.getString(),
			    cacheKey, bot, sourceMeshRequest);
		    if (assetPath)
			bu_free(assetPath, "stream LoD asset path");
		} else {
		    generatedGeometry = cad_mesh_part_geometry_from_internal(
			localIntern, data->source, generated);
		}
		if (generatedGeometry &&
		    drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE)
		    (void)cad_mesh_append_hidden_line_edges(generated);
		if (generatedGeometry) {
		    const char *typeLabel = primitive_type_label(localIntern);
		    if (!transformedReuse)
			data->cache->storeMeshCadGeometry(cacheKey,
			    std::move(generated), typeLabel,
			    lodBacked ? "aabb" : "surface", NULL,
			    lodBacked,
			    lodBacked ? &sourceMeshRequest : NULL);
		    cachedMesh = data->cache->findMeshCadGeometry(cacheKey);
		}
	    }
	}
	if (cachedWire) {
	    bobol_performance_counter_add(hadCachedWire ?
		BOBOL_PERF_MESH_CACHE_HITS : BOBOL_PERF_MESH_CACHE_MISSES,
		1);
	    char *path = db_path_to_string(pathp);
	    const char *typeLabel = cachedWire->sourceType.empty() ? "wire" :
		cachedWire->sourceType.c_str();
	    const char *geometryKind = cachedWire->geometryKind.empty() ? "line" :
		cachedWire->geometryKind.c_str();
	    compact_occurrence_build input;
	    input.occurrence.geometry = cachedWire->geometry;
	    input.occurrence.viewDependentCsgGeometry =
		cachedWire->viewDependentCsgGeometry ? TRUE : FALSE;
	    input.occurrence.localTransform = mat_to_sbmatrix(tsp->ts_mat);
	    input.occurrence.summary = compact_occurrence_tree_summary(
		data->source, tsp, pathp, path, dp->d_namep,
		geometryKind && BU_STR_EQUAL(geometryKind, "annotation") ?
		"annotation" : typeLabel,
		geometryKind, data->revision,
		BObolRealizedShapeSummary::SHAPE_VLIST,
		static_cast<BObolMaterialColorSweep *>(data->material_sweep));
	    input.occurrence.occurrenceIndex =
		data->source->occurrenceIndex.getValue();
	    input.occurrence.booleanOperation =
		data->source->booleanOperation.getValue();
	    input.semantic = compact_semantic_from_summary(
		input.occurrence.summary);
	    input.dashed = (tsp->ts_sofar & TS_SOFAR_MINUS) ? TRUE : FALSE;
	    const size_t entryCount = data->compact_index->entries.size();
	    compact_add_occurrence(data->source, *data->compact_index, input,
		data->compact_ordinal, data->compact_unsupported);
	    compact_apply_walk_identity(data->source, *data->compact_index,
		entryCount, tsp, pathp, walkOccurrenceIdentity, duplicateOrdinal);
	    if (data->compact_index->entries.size() > entryCount)
		realize_walk_stream_push_current(data, input.occurrence,
		    *data->compact_index, entryCount);
	    SbBox3f bounds = database_source_transform_bounds(
		compact_part_geometry_bounds(cachedWire->geometry),
		input.occurrence.localTransform);
	    realize_walk_extend_bounds(data, bounds);
	    data->realized_shapes++;
	    if (path)
		bu_free(path, "db_path_to_string");
	    return make_nop_tree();
	}
    }

    SoBRLVListShape *sharedVListShape = NULL;
    SoBRLMeshShape *sharedMeshShape = NULL;
    SbBox3f cacheBounds;
    std::string cacheKey = realize_geometry_cache_key(dp);
    source_lod_cache_key_append(cacheKey, data->source, cacheBounds,
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP);
    BObolRealizationCacheMap<SoBRLVListShape *>::iterator foundVList =
	data->cache->sharedMeshVListGeometry.find(cacheKey);
    if (foundVList != data->cache->sharedMeshVListGeometry.end()) {
	sharedVListShape = foundVList->second;
    } else {
	BObolRealizationCacheMap<SoBRLMeshShape *>::iterator foundMesh =
	    data->cache->sharedMeshGeometry.find(cacheKey);
	if (foundMesh != data->cache->sharedMeshGeometry.end())
	    sharedMeshShape = foundMesh->second;
    }
    const BObolCachedPartGeometry *cachedWire = data->compact_index ?
	data->cache->findMeshVListCadGeometry(cacheKey) : NULL;
    const BObolCachedPartGeometry *cachedMesh = data->compact_index ?
	data->cache->findMeshCadGeometry(cacheKey) : NULL;
    if (sharedVListShape &&
	!source_cached_wire_matches_mesh_presentation(data->source,
	    sharedVListShape->sourceType.getValue().getString(),
	    sharedVListShape->geometryKind.getValue().getString()))
	sharedVListShape = NULL;
    if (cachedWire && !source_cached_wire_matches_mesh_presentation(
	    data->source, cachedWire->sourceType.c_str(),
	    cachedWire->geometryKind.c_str()))
	cachedWire = NULL;
    if (!source_cached_mesh_matches_presentation(data->source, dp)) {
	sharedMeshShape = NULL;
	cachedMesh = NULL;
    }
    bobol_performance_counter_add(
	(sharedVListShape || sharedMeshShape || cachedWire || cachedMesh) ?
	BOBOL_PERF_MESH_CACHE_HITS :
	BOBOL_PERF_MESH_CACHE_MISSES, 1);

    const char *typeLabel = sharedVListShape ?
			    sharedVListShape->sourceType.getValue().getString() :
			    (sharedMeshShape ?
			     sharedMeshShape->sourceType.getValue().getString() :
			     (cachedWire && !cachedWire->sourceType.empty() ?
			      cachedWire->sourceType.c_str() :
			      (cachedMesh && !cachedMesh->sourceType.empty() ?
			       cachedMesh->sourceType.c_str() : NULL)));
    if (!sharedVListShape && !sharedMeshShape && !cachedWire && !cachedMesh) {
	owned_leaf_internal validInternal;
	struct rt_db_internal *localIntern =
	    import_walk_leaf_internal(tsp, dp, &validInternal);
	if (!localIntern) {
	    data->failed_shapes++;
	    set_leaf_import_diagnostic(data, pathp,
		validInternal.ownsLocal ? &validInternal.local : NULL);
	    return TREE_NULL;
	}
	if (database_source_terminal_empty_bot(localIntern))
	    return make_nop_tree();

	typeLabel = primitive_type_label(localIntern);
	const int internalType = localIntern->idb_type;
	SbBox3f localBounds;
	if (source_view_lod_active(data->source))
	    (void)local_bounds_from_internal(localIntern, localBounds);
	if (internalType == ID_MATERIAL) {
	    SoBRLMaterialObject *materialObject =
		material_object_from_internal(static_cast<struct rt_material_internal *>(localIntern->idb_ptr));
	    if (!materialObject) {
		data->failed_shapes++;
		set_walk_diagnostic(data, pathp,
				    "material object realization failed");
		return TREE_NULL;
	    }
	    char *path = db_path_to_string(pathp);
	    SoSeparator *leaf = new SoSeparator;
	    assign_material_identity(materialObject, path, dp->d_namep,
				     typeLabel, data->revision);
	    leaf->addChild(materialObject);
	    database_source_add_realized_child(data->source, leaf);
	    data->realized_shapes++;
	    if (path)
		bu_free(path, "db_path_to_string");
	    return make_nop_tree();
	}

	if (internalType == ID_PNTS) {
	    sharedVListShape = vlist_from_pnts(
				   static_cast<const struct rt_pnts_internal *>(localIntern->idb_ptr));
	    if (sharedVListShape)
		assign_shared_geometry_identity(sharedVListShape,
						dp->d_namep, typeLabel, data->revision, "point");
	} else if (primitive_uses_wire_in_mesh_mode(internalType)) {
	    sharedVListShape = vlist_from_lod_realization_internal(localIntern,
		data->source, localBounds);
	    if (!sharedVListShape)
		sharedVListShape = vlist_from_plot_internal(localIntern,
		    data->source);
	    if (sharedVListShape) {
		assign_shared_geometry_identity(sharedVListShape,
						dp->d_namep, typeLabel, data->revision, "line");
		if (primitive_is_annotation(internalType, typeLabel)) {
		    sharedVListShape->sourceType = "annotation";
		    sharedVListShape->geometryKind = "annotation";
		}
	    }
	} else if ((source_record_draw_mode(data->source) ==
		    BOBOL_LOD_DRAW_WIRE ||
		    (source_record_draw_mode(data->source) ==
			BOBOL_LOD_DRAW_SHADED_BOTS &&
		     (!localIntern->idb_meth ||
		      !localIntern->idb_meth->ft_indexed_face_set))) &&
			   internalType != ID_BOT) {
	    sharedVListShape = vlist_from_lod_realization_internal(localIntern,
		data->source, localBounds);
	    if (!sharedVListShape)
		sharedVListShape = vlist_from_plot_internal(localIntern,
		    data->source);
	    if (sharedVListShape)
		assign_shared_geometry_identity(sharedVListShape,
						dp->d_namep, typeLabel, data->revision, "line");
	} else {
	    /* Non-BOT mesh primitives realized on a cache miss.  BOTs are
	     * converted directly to compact geometry by the (parallel) prefill
	     * and hit the cache above, so they do not reach this build path. */
	    sharedMeshShape = mesh_from_internal(localIntern, data->source);
	    if (sharedMeshShape)
		assign_shared_geometry_identity(sharedMeshShape,
				dp->d_namep, typeLabel, data->revision, "surface");
	}

	if (!sharedVListShape && !sharedMeshShape) {
	    char reason[256] = {0};
	    data->failed_shapes++;
	    snprintf(reason, sizeof(reason),
		     "unsupported or failed mesh conversion/tessellation for primitive type '%s'",
		     typeLabel);
	    set_walk_diagnostic(data, pathp, reason);
	    return TREE_NULL;
	}

	if (sharedVListShape) {
	    if (!data->compact_index)
		data->cache->storeMeshVListGeometry(cacheKey, sharedVListShape);
	    typeLabel = sharedVListShape->sourceType.getValue().getString();
	} else {
	    if (!data->compact_index)
		data->cache->storeMeshGeometry(cacheKey, sharedMeshShape);
	    typeLabel = sharedMeshShape->sourceType.getValue().getString();
	}
    }

    char *path = db_path_to_string(pathp);
    if (data->compact_index) {
	const SbMatrix localMatrix = mat_to_sbmatrix(tsp->ts_mat);
	const size_t entryCount = data->compact_index->entries.size();
	compact_occurrence_build input;
	input.occurrence.localTransform = localMatrix;
	input.occurrence.occurrenceIndex =
	    data->source->occurrenceIndex.getValue();
	input.occurrence.booleanOperation =
	    data->source->booleanOperation.getValue();
	if (sharedVListShape || cachedWire) {
	    const char *geometryKind = sharedVListShape ?
		sharedVListShape->geometryKind.getValue().getString() :
		(cachedWire && !cachedWire->geometryKind.empty() ?
		 cachedWire->geometryKind.c_str() : "line");
	    input.occurrence.geometry =
		cachedWire ? cachedWire->geometry :
		std::shared_ptr<const Obol::PartGeometry>();
	    input.occurrence.viewDependentCsgGeometry = cachedWire &&
		cachedWire->viewDependentCsgGeometry ? TRUE : FALSE;
	    if (!input.occurrence.geometry && sharedVListShape) {
		Obol::PartGeometryBuilder generated;
		if (cad_vlist_part_geometry(sharedVListShape, generated))
		    input.occurrence.geometry =
			data->cache->storeMeshVListCadGeometry(cacheKey,
			    std::move(generated), typeLabel, geometryKind);
	    }
	    input.occurrence.summary = compact_occurrence_tree_summary(
		data->source, tsp, pathp, path, dp->d_namep,
		geometryKind && BU_STR_EQUAL(geometryKind, "annotation") ?
		"annotation" : typeLabel,
		geometryKind, data->revision,
		BObolRealizedShapeSummary::SHAPE_VLIST,
		static_cast<BObolMaterialColorSweep *>(data->material_sweep));
	} else {
	    const char *geometryKind = sharedMeshShape ?
		sharedMeshShape->geometryKind.getValue().getString() :
		(cachedMesh && !cachedMesh->geometryKind.empty() ?
		 cachedMesh->geometryKind.c_str() : "surface");
	    if (sharedMeshShape) {
		sharedMeshShape->drawMode = source_record_draw_mode(data->source);
		sharedMeshShape->hiddenLine =
		    source_record_draw_mode(data->source) ==
		    BOBOL_LOD_DRAW_HIDDEN_LINE ? TRUE : FALSE;
		if (typeLabel && BU_STR_EQUAL(typeLabel, "bot") &&
		    data->source->lodBotThreshold.getValue() > 0)
		    publish_lod_mesh_if_available(sharedMeshShape, data->source,
			tsp->ts_dbip, dp->d_namep);
	    }
	    input.occurrence.geometry =
		cachedMesh ? cachedMesh->geometry :
		std::shared_ptr<const Obol::PartGeometry>();
	    if (cachedMesh)
		input.occurrence.geometryTransform = cachedMesh->geometryTransform;
	    if (!input.occurrence.geometry && sharedMeshShape) {
		Obol::PartGeometryBuilder generated;
		if (cad_mesh_part_geometry(sharedMeshShape, generated))
		    input.occurrence.geometry =
			data->cache->storeMeshCadGeometry(cacheKey,
			    std::move(generated), typeLabel, geometryKind, NULL,
			    sharedMeshShape->isLodBackedMesh());
	    }
	    input.occurrence.lodBacked = sharedMeshShape ?
		sharedMeshShape->isLodBackedMesh() :
		(cachedMesh && cachedMesh->lodBacked ? TRUE : FALSE);
	    if (sharedMeshShape)
		input.occurrence.sourceMeshRequestValid =
		    sharedMeshShape->makeSourceMeshRequest(
			input.occurrence.sourceMeshRequest);
	    else if (cachedMesh && cachedMesh->sourceMeshRequestValid) {
		input.occurrence.sourceMeshRequestValid = TRUE;
		input.occurrence.sourceMeshRequest =
		    cachedMesh->sourceMeshRequest;
	    }
	    input.occurrence.summary = compact_occurrence_tree_summary(
		data->source, tsp, pathp, path, dp->d_namep, typeLabel,
		geometryKind,
		data->revision, BObolRealizedShapeSummary::SHAPE_MESH,
		static_cast<BObolMaterialColorSweep *>(data->material_sweep));
	    if (sharedMeshShape) {
		BObolRealizedShapeSummary meshSummary;
		realized_mesh_shape_summary(sharedMeshShape, meshSummary);
		input.occurrence.summary.lodAvailable = meshSummary.lodAvailable;
		input.occurrence.summary.lodActiveCut =
		    meshSummary.lodActiveCut;
		input.occurrence.summary.lodFaceCount = meshSummary.lodFaceCount;
		input.occurrence.summary.lodPointCount = meshSummary.lodPointCount;
		input.occurrence.summary.lodOriginalPointCount =
		    meshSummary.lodOriginalPointCount;
		input.occurrence.summary.lodNormalCount =
		    meshSummary.lodNormalCount;
		input.occurrence.summary.lodHasSnappedPoints =
		    meshSummary.lodHasSnappedPoints;
		input.occurrence.summary.lodHasNormals = meshSummary.lodHasNormals;
		input.occurrence.summary.lodBoundsMin = meshSummary.lodBoundsMin;
		input.occurrence.summary.lodBoundsMax = meshSummary.lodBoundsMax;
	    }
	    if (input.occurrence.sourceMeshRequestValid)
		compact_source_mesh_request_sync(
		    input.occurrence.sourceMeshRequest,
		    input.occurrence.summary);
	    if (input.occurrence.sourceMeshRequestValid)
		compact_summary_lod_from_source_mesh_request(
		    input.occurrence.summary, input.occurrence.sourceMeshRequest);
	    if (input.occurrence.sourceMeshRequestValid)
		cache_mesh_cad_source_request(data->cache, cacheKey,
		    input.occurrence.sourceMeshRequest);
	}
	if (!input.occurrence.geometry) {
	    data->compact_unsupported = 1;
	} else {
	    input.semantic = compact_semantic_from_summary(
		input.occurrence.summary);
	    compact_add_occurrence(data->source, *data->compact_index, input,
		data->compact_ordinal, data->compact_unsupported);
	    compact_apply_walk_identity(data->source, *data->compact_index,
		entryCount, tsp, pathp, walkOccurrenceIdentity,
		duplicateOrdinal);
	    if (data->compact_index->entries.size() > entryCount)
		realize_walk_stream_push_current(data, input.occurrence,
		    *data->compact_index, entryCount);
	}
	BObolCompactInstanceEntry *entry =
	    data->compact_index->entries.size() > entryCount ?
	    &data->compact_index->entries.back() : NULL;
	SbBox3f bounds = entry ?
	    database_source_transform_bounds(
		compact_part_geometry_bounds(entry->geometry),
		entry->localTransform) :
	    SbBox3f();
	realize_walk_extend_bounds(data, bounds);
    } else {
	SoSeparator *leaf = realize_instance_leaf_separator(tsp);
	if (sharedVListShape) {
	    SoBRLVListShape *vlistShape = new SoBRLVListShape;
	    assign_realized_identity(vlistShape, tsp, path, dp->d_namep,
		typeLabel, data->revision, data->source,
	    static_cast<BObolMaterialColorSweep *>(data->material_sweep));
	    vlistShape->setSharedGeometry(sharedVListShape);
	    const char *geometryKind =
		sharedVListShape->geometryKind.getValue().getString();
	    vlistShape->geometryKind = geometryKind && geometryKind[0] ?
		geometryKind : "line";
	    if (geometryKind && BU_STR_EQUAL(geometryKind, "annotation"))
		vlistShape->sourceType = "annotation";
	    leaf->addChild(vlistShape);
	} else {
	    SoBRLMeshShape *shape = mesh_instance_for_shared_geometry(sharedMeshShape);
	    assign_realized_identity(shape, tsp, path, dp->d_namep, typeLabel,
		data->revision, data->source,
	    static_cast<BObolMaterialColorSweep *>(data->material_sweep));
	    shape->setSharedGeometry(sharedMeshShape);
	    const char *geometryKind =
		sharedMeshShape->geometryKind.getValue().getString();
	    shape->geometryKind = geometryKind && geometryKind[0] ?
		geometryKind : "surface";
	    if (typeLabel && BU_STR_EQUAL(typeLabel, "bot") &&
		data->source->lodBotThreshold.getValue() > 0)
		publish_lod_mesh_if_available(shape, data->source, tsp->ts_dbip,
		    dp->d_namep);
	    leaf->addChild(shape);
	}
	database_source_add_realized_child(data->source, leaf);
    }
    data->realized_shapes++;
    if (path)
	bu_free(path, "db_path_to_string");

    return make_nop_tree();
}


static void
source_bounds_extend_points(SbBox3f &bounds,
			    const SoMFVec3f &points,
			    const SbMatrix &matrix)
{
    for (int i = 0; i < points.getNum(); i++) {
	SbVec3f transformed;
	matrix.multVecMatrix(points[i], transformed);
	bounds.extendBy(transformed);
    }
}


static SbBool
source_bounds_for_realized_node(const SoNode *node,
				const SbMatrix &matrix,
				SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!node || node_is_source_placement_transform(node) ||
	node_is_auxiliary_vlist(node) || node_is_auxiliary_source(node))
	return FALSE;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	const SoBRLVListShape *shape =
	    static_cast<const SoBRLVListShape *>(node);
	const SoBRLVListShape *geom = shape->getGeometrySource();
	if (!geom || geom->point.getNum() <= 0)
	    return FALSE;
	source_bounds_extend_points(bounds, geom->point, matrix);
	return bounds.isEmpty() ? FALSE : TRUE;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	const SoBRLMeshShape *shape =
	    static_cast<const SoBRLMeshShape *>(node);
	const SoBRLMeshShape *geom = shape->getGeometrySource();
	if (!geom || geom->point.getNum() <= 0)
	    return FALSE;
	source_bounds_extend_points(bounds, geom->point, matrix);
	return bounds.isEmpty() ? FALSE : TRUE;
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return FALSE;

    const SoGroup *group = static_cast<const SoGroup *>(node);
    SbMatrix childMatrix = matrix;
    SbBool valid = FALSE;
    for (int i = 0; i < group->getNumChildren(); i++) {
	const SoNode *child = group->getChild(i);
	if (!child)
	    continue;
	if (child->isOfType(SoMatrixTransform::getClassTypeId())) {
	    const SoMatrixTransform *transform =
		static_cast<const SoMatrixTransform *>(child);
	    childMatrix.multRight(transform->matrix.getValue());
	    continue;
	}

	SbBox3f childBounds;
	if (source_bounds_for_realized_node(child, childMatrix,
					    childBounds)) {
	    bounds.extendBy(childBounds);
	    valid = TRUE;
	}
    }

    return valid && !bounds.isEmpty() ? TRUE : FALSE;
}


static void
update_source_bounds_from_realized_geometry(SoBRLDatabaseSource *source,
	SbBool exact = FALSE)
{
    if (!source)
	return;
    if (!exact && source->hasExactSourceBounds())
	return;

    SbBox3f bounds;
    bounds.makeEmpty();
    SbBool valid = FALSE;
    const SbMatrix identity = SbMatrix::identity();
    for (int i = 0; i < source->getNumChildren(); i++) {
	SbBox3f childBounds;
	if (source_bounds_for_realized_node(source->getChild(i), identity,
					    childBounds)) {
	    bounds.extendBy(childBounds);
	    valid = TRUE;
	}
    }

    if (valid && !bounds.isEmpty()) {
	(void)source->setSourceBoundsState(TRUE, bounds.getMin(),
					   bounds.getMax(), exact);
    } else {
	source->clearSourceBounds();
    }
}

static int
source_has_auxiliary_children(const SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;

    for (int i = 0; i < source->getNumChildren(); i++) {
	SoNode *child = source->getChild(i);
	if (node_is_auxiliary_vlist(child) || node_is_auxiliary_source(child))
	    return 1;
    }
    return 0;
}

static void
cad_shape_color(const SbBool selected,
		const SbColor &selectedColor,
		const SbBool highlighted,
		const SbColor &highlightedColor,
		const SbBool ghosted,
		const SbColor &ghostedColor,
		const SbBool colorOverride,
		const SbColor &overrideColor,
		const SbBool materialColorValid,
		const SbColor &materialColor,
		const SbColor &fallbackColor,
		float transparency,
		SbColor4f &colorOut)
{
    SbColor color = fallbackColor;
    float alpha = 1.0f;

    if (highlighted) {
	color = highlightedColor;
    } else if (selected) {
	color = selectedColor;
    } else if (ghosted) {
	color = ghostedColor;
	alpha = 0.35f;
    } else if (colorOverride) {
	color = overrideColor;
    } else if (materialColorValid) {
	color = materialColor;
    }

    if (transparency > 0.0f) {
	if (transparency > 1.0f)
	    transparency = 1.0f;
	alpha *= (1.0f - transparency);
    }

    colorOut = SbColor4f(color[0], color[1], color[2], alpha);
}

uint64_t
compact_next_revision(uint64_t revision)
{
    return bobol_identity_successor_or_terminate(revision);
}

bool
compact_style_equal(const Obol::InstanceStyle &a,
	const Obol::InstanceStyle &b)
{
    return a.hasColorOverride == b.hasColorOverride &&
	!database_source_float_different(a.color[0], b.color[0]) &&
	!database_source_float_different(a.color[1], b.color[1]) &&
	!database_source_float_different(a.color[2], b.color[2]) &&
	!database_source_float_different(a.color[3], b.color[3]) &&
	!database_source_float_different(a.lineWidth, b.lineWidth) &&
	a.linePattern == b.linePattern &&
	a.linePatternFactor == b.linePatternFactor;
}

bool
compact_semantic_equal(const SoBRLCadAssembly::InstanceSemantic &a,
    const SoBRLCadAssembly::InstanceSemantic &b)
{
    return a.path == b.path &&
	a.sourceInstanceKey == b.sourceInstanceKey &&
	a.sourceName == b.sourceName &&
	a.sourceType == b.sourceType &&
	a.materialShader == b.materialShader &&
	a.editIntentId == b.editIntentId &&
	a.editIntentRole == b.editIntentRole &&
	a.sourceId == b.sourceId && a.regionId == b.regionId &&
	a.airCode == b.airCode && a.materialId == b.materialId &&
	a.los == b.los &&
	a.materialColorValid == b.materialColorValid &&
	(!a.materialColorValid ||
	 database_source_color_equal(a.materialColor, b.materialColor)) &&
	a.primitiveKind == b.primitiveKind;
}

void
compact_note_semantic_change(BObolCompactInstanceEntry &entry)
{
    entry.semanticRevision = compact_next_revision(entry.semanticRevision);
}

Obol::InstanceStyle
compact_effective_style(const BObolCompactInstanceEntry &entry)
{
    return compact_presentation_style(entry);
}

bool
compact_retired_overview(const BObolCompactInstanceEntry &entry)
{
    return !entry.authoredVisible &&
	BU_STR_EQUAL(entry.shapeSummary.recordRole.getString(), "lod-overview");
}

SbBool
compact_effective_authored_visibility(
    const BObolCompactInstanceEntry &entry)
{
    return compact_presentation_visibility(entry, compact_retired_overview(entry));
}

SbBool
compact_effective_highlight(const BObolCompactInstanceEntry &entry)
{
    return compact_presentation_highlight(entry);
}

void
compact_sync_shape_display_summary(BObolCompactInstanceEntry &entry)
{
    auto &summary = entry.shapeSummary;
    summary.visible = entry.visible;
    summary.selectable = entry.selectable;
    summary.selected = entry.selected;
    summary.highlighted = entry.highlighted;
    summary.lineStyle = entry.style.linePattern == 0xffffu ? 0 : 1;
    summary.lineWidth = entry.style.lineWidth > 0.0f ?
	static_cast<int>(entry.style.lineWidth + 0.5f) : 0;
    summary.transparency = std::max(0.0f,
	std::min(1.0f, 1.0f - entry.style.color[3]));
}

void
compact_sync_shape_summary_state(BObolCompactInstanceEntry &entry)
{
    BObolRealizedShapeSummary &summary = entry.shapeSummary;
    summary.valid = TRUE;
    /* shapeKind describes the source primitive's semantic role, while the
     * geometry flags describe the channels that can be drawn right now.  A
     * progressive mesh initially owns only a wire proxy, but it must remain a
     * mesh for LoD submission, picking identity, and exact export. */
    summary.shapeKind = (entry.meshGeometry ||
	(entry.lodBacked && entry.sourceMeshRequestValid)) ?
	BObolRealizedShapeSummary::SHAPE_MESH :
	BObolRealizedShapeSummary::SHAPE_VLIST;
    summary.path = entry.semantic.path;
    summary.sourceName = entry.semantic.sourceName;
    summary.sourceType = entry.semantic.sourceType;
    summary.sourceId = entry.semantic.sourceId;
    summary.regionId = entry.semantic.regionId;
    summary.airCode = entry.semantic.airCode;
    summary.materialId = entry.semantic.materialId;
    summary.los = entry.semantic.los;
    summary.materialColorValid = entry.semantic.materialColorValid;
    summary.materialColor = entry.semantic.materialColor;
    summary.materialShader = entry.semantic.materialShader;
    summary.editIntentId = entry.semantic.editIntentId;
    summary.editIntentRole = entry.semantic.editIntentRole;
    compact_sync_shape_display_summary(entry);
}

static void
compact_sync_shape_summary(BObolCompactInstanceEntry &entry)
{
    compact_sync_shape_summary_state(entry);
    BObolRealizedShapeSummary &summary = entry.shapeSummary;
    summary.bounds = database_source_transform_bounds(
	compact_part_geometry_bounds(entry.geometry), entry.geometryTransform);
    summary.boundsValid = !summary.bounds.isEmpty();
    summary.pointCount = 0;
    summary.commandCount = 0;
    summary.segmentCount = 0;
    summary.triangleCount = 0;
    summary.indexCount = 0;
    if (entry.geometry && entry.geometry->points) {
	const Obol::PointRep &points = *entry.geometry->points;
	summary.pointCount = static_cast<int>(points.positions.size());
	summary.commandCount = summary.pointCount;
	summary.pointPrimitiveCount = summary.pointCount;
    }
    if (entry.geometry && entry.geometry->wire) {
	const Obol::WireRep &wire = *entry.geometry->wire;
	const int wirePointCount = static_cast<int>(wire.segmentPoints.size());
	summary.pointCount += wirePointCount;
	summary.commandCount += wirePointCount;
	summary.segmentCount = static_cast<int>(wire.segmentCount());
	for (const Obol::WirePolyline &polyline : wire.polylines) {
	    summary.pointCount += static_cast<int>(polyline.points.size());
	    summary.segmentCount += polyline.points.empty() ? 0 :
		static_cast<int>(polyline.points.size() - 1);
	}
    }
    if (entry.geometry && entry.geometry->shaded) {
	const Obol::TriMesh &mesh = *entry.geometry->shaded;
	summary.pointCount = static_cast<int>(mesh.positions.size());
	summary.indexCount = static_cast<int>(mesh.indices.size());
	summary.triangleCount = static_cast<int>(mesh.indices.size() / 3);
    }
}

static bool
compact_appearance_equal(const BObolCompactInstanceEntry &a,
	const BObolCompactInstanceEntry &b)
{
    return compact_style_equal(a.normalStyle, b.normalStyle) &&
	compact_style_equal(a.selectedStyle, b.selectedStyle) &&
	compact_style_equal(a.highlightedStyle, b.highlightedStyle);
}

static Obol::InstanceStyle
cad_vlist_style_state(const SoBRLVListShape *shape, SbBool selected,
	SbBool highlighted)
{
    Obol::InstanceStyle style;
    if (!shape)
	return style;

    style.hasColorOverride = true;
    cad_shape_color(selected, shape->selectedColor.getValue(), highlighted,
		    shape->highlightedColor.getValue(),
		    shape->ghosted.getValue(), shape->ghostedColor.getValue(),
		    shape->colorOverride.getValue(), shape->color.getValue(),
		    shape->materialColorValid.getValue(),
		    shape->materialColor.getValue(), shape->color.getValue(),
		    shape->transparency.getValue(), style.color);
    style.lineWidth = shape->lineWidth.getValue() > 0 ?
		      static_cast<float>(shape->lineWidth.getValue()) : 1.0f;
    if (shape->lineStyle.getValue() != 0)
	style.linePattern = 0xcf33u;
    return style;
}

static Obol::InstanceStyle
cad_vlist_style(const SoBRLVListShape *shape)
{
    return cad_vlist_style_state(shape,
	shape ? shape->selected.getValue() : FALSE,
	shape ? shape->highlighted.getValue() : FALSE);
}

static Obol::InstanceStyle
cad_mesh_style_state(const SoBRLMeshShape *shape, SbBool selected,
	SbBool highlighted)
{
    Obol::InstanceStyle style;
    if (!shape)
	return style;

    style.hasColorOverride = true;
    cad_shape_color(selected, shape->selectedColor.getValue(), highlighted,
		    shape->highlightedColor.getValue(),
		    shape->ghosted.getValue(), shape->ghostedColor.getValue(),
		    shape->colorOverride.getValue(), shape->color.getValue(),
		    shape->materialColorValid.getValue(),
		    shape->materialColor.getValue(), shape->color.getValue(),
		    shape->transparency.getValue(), style.color);
    style.lineWidth = shape->lineWidth.getValue() > 0 ?
		      static_cast<float>(shape->lineWidth.getValue()) : 1.0f;
    if (shape->lineStyle.getValue() != 0)
	style.linePattern = 0xcf33u;
    return style;
}

static Obol::InstanceStyle
cad_mesh_style(const SoBRLMeshShape *shape)
{
    return cad_mesh_style_state(shape,
	shape ? shape->selected.getValue() : FALSE,
	shape ? shape->highlighted.getValue() : FALSE);
}

static int
cad_vlist_point(const SoBRLVListShape *shape,
		const SoBRLVListShape *geom,
		int index,
		SbVec3f &point)
{
    if (!geom || index < 0 || index >= geom->point.getNum())
	return 0;

    double precise[3] = {0.0, 0.0, 0.0};
    if (shape && shape->getPrecisePoint(index, precise)) {
	point.setValue(static_cast<float>(precise[0]),
		       static_cast<float>(precise[1]),
		       static_cast<float>(precise[2]));
	return 1;
    }

    point = geom->point[index];
    return 1;
}

class CadGeometryHash {
public:
    CadGeometryHash(void) :
	first(14695981039346656037ULL),
	second(7809847782465536322ULL)
    {
    }

    void appendByte(uint8_t value)
    {
	first = (first ^ value) * 1099511628211ULL;
	second = (second ^ value) * 14029467366897019727ULL;
    }

    void appendU32(uint32_t value)
    {
	first = (first ^ value) * 1099511628211ULL;
	second = (second ^ value) * 14029467366897019727ULL;
    }

    void appendU64(uint64_t value)
    {
	this->appendU32(static_cast<uint32_t>(value));
	this->appendU32(static_cast<uint32_t>(value >> 32));
    }

    void appendFloat(float value)
    {
	uint32_t bits = 0;
	static_assert(sizeof(bits) == sizeof(value), "float hash encoding");
	memcpy(&bits, &value, sizeof(bits));
	this->appendU32(bits);
    }

    void appendVec3(const SbVec3f &value)
    {
	this->appendFloat(value[0]);
	this->appendFloat(value[1]);
	this->appendFloat(value[2]);
    }

    void appendBox(const SbBox3f &value)
    {
	this->appendByte(value.isEmpty() ? 0 : 1);
	if (value.isEmpty())
	    return;
	this->appendVec3(value.getMin());
	this->appendVec3(value.getMax());
    }

    void appendQuantization(const Obol::ProgressiveQuantization &value)
    {
	this->appendByte(value.xBits);
	this->appendByte(value.yBits);
	this->appendByte(value.zBits);
    }

    void appendString(const char *value)
    {
	if (!value)
	    return;
	for (; *value; value++)
	    this->appendByte(static_cast<uint8_t>(*value));
	this->appendByte(0);
    }

    Obol::PartId id(void) const
    {
	Obol::PartId result;
	result.w0 = first;
	result.w1 = second;
	return result;
    }

private:
    uint64_t first;
    uint64_t second;
};

template <typename Geometry>
static int
cad_part_key_for_geometry(const char *kind,
			  const Geometry &geometry,
			  std::string &key)
{
    if (!kind)
	return 0;

    CadGeometryHash hash;
    hash.appendString(kind);
    /* Part identity covers presentation semantics as well as vertex arrays.
     * A degenerate AABB around a line can have byte-identical endpoints to
     * the authored line it temporarily represents.  Deduplicating those two
     * records under one PartId leaves whichever structuralProxy marker was
     * inserted first authoritative for every occurrence, so a completed
     * stream may continue drawing and counting the box forever. */
    hash.appendByte(geometry.shadedCullBackfaces ? 1 : 0);
    hash.appendByte(geometry.subpixelProxyEligible ? 1 : 0);
    hash.appendByte(geometry.structuralProxy ? 1 : 0);
    hash.appendByte(geometry.conservativeBounds ? 1 : 0);
    if (geometry.conservativeBounds)
	hash.appendBox(*geometry.conservativeBounds);
    hash.appendByte(geometry.points ? 1 : 0);
    if (geometry.points) {
	const Obol::PointRep &points = *geometry.points;
	hash.appendU32(
	    static_cast<uint32_t>(points.positions.size()));
	for (const SbVec3f &point : points.positions)
	    hash.appendVec3(point);
	hash.appendU32(
	    static_cast<uint32_t>(points.pointIds.size()));
	for (uint32_t id : points.pointIds)
	    hash.appendU32(id);
	hash.appendU32(
	    static_cast<uint32_t>(points.colorValid.size()));
	for (uint8_t valid : points.colorValid)
	    hash.appendByte(valid);
	hash.appendU32(
	    static_cast<uint32_t>(points.colors.size()));
	for (const SbColor &color : points.colors)
	    hash.appendVec3(color);
	hash.appendU32(
	    static_cast<uint32_t>(points.scaleValid.size()));
	for (uint8_t valid : points.scaleValid)
	    hash.appendByte(valid);
	hash.appendU32(
	    static_cast<uint32_t>(points.scales.size()));
	for (float scale : points.scales)
	    hash.appendFloat(scale);
	hash.appendU32(
	    static_cast<uint32_t>(points.normalValid.size()));
	for (uint8_t valid : points.normalValid)
	    hash.appendByte(valid);
	hash.appendU32(
	    static_cast<uint32_t>(points.normals.size()));
	for (const SbVec3f &normal : points.normals)
	    hash.appendVec3(normal);
	hash.appendBox(points.bounds);
    }
    hash.appendByte(geometry.wire ? 1 : 0);
    if (geometry.wire) {
	const Obol::WireRep &wire = *geometry.wire;
	hash.appendU32(
	    static_cast<uint32_t>(wire.segmentPoints.size()));
	for (const SbVec3f &point : wire.segmentPoints)
	    hash.appendVec3(point);
	hash.appendU32(
	    static_cast<uint32_t>(wire.segmentIds.size()));
	for (uint32_t id : wire.segmentIds)
	    hash.appendU32(id);
	hash.appendU32(
	    static_cast<uint32_t>(wire.polylines.size()));
	for (const Obol::WirePolyline &polyline : wire.polylines) {
	    hash.appendU32(polyline.edgeId);
	    hash.appendU32(
		static_cast<uint32_t>(polyline.points.size()));
	    for (const SbVec3f &point : polyline.points)
		hash.appendVec3(point);
	}
	hash.appendBox(wire.bounds);
	const bool progressiveWire = wire.isProgressive();
	hash.appendByte(progressiveWire ? 1 : 0);
	if (progressiveWire) {
	    hash.appendU32(static_cast<uint32_t>(wire.progressiveCuts.size()));
	    for (const Obol::ProgressiveWireCut &cut : wire.progressiveCuts) {
		hash.appendU32(cut.segmentFirst);
		hash.appendU32(cut.segmentCount);
		hash.appendQuantization(cut.quantization);
		hash.appendFloat(cut.maximumNormalizedError);
	    }
	    hash.appendByte(wire.progressiveMinimumCut);
	    hash.appendByte(wire.progressiveResidentCut);
	    hash.appendVec3(wire.progressiveQuantizationMinimum);
	    hash.appendVec3(wire.progressiveQuantizationMaximum);
	    hash.appendU64(wire.progressiveLineage);
	    hash.appendU32(static_cast<uint32_t>(
		wire.progressiveClusters.size()));
	    hash.appendU32(wire.progressiveClusterGridResolution);
	    for (const Obol::ProgressiveWireCluster &cluster :
		 wire.progressiveClusters) {
		hash.appendBox(cluster.bounds);
		hash.appendByte(cluster.residentCut);
		hash.appendU32(static_cast<uint32_t>(cluster.ranges.size()));
		for (const Obol::ProgressiveWireClusterRange &range :
		     cluster.ranges) {
		    hash.appendU32(range.firstSegment);
		    hash.appendU32(range.segmentCount);
		    hash.appendByte(range.activationCut);
		}
	    }
	}
    }
    hash.appendByte(geometry.shaded ? 1 : 0);
    if (geometry.shaded) {
	const Obol::TriMesh &mesh = *geometry.shaded;
	hash.appendU32(
	    static_cast<uint32_t>(mesh.positions.size()));
	for (const SbVec3f &point : mesh.positions)
	    hash.appendVec3(point);
	hash.appendU32(
	    static_cast<uint32_t>(mesh.normals.size()));
	for (const SbVec3f &normal : mesh.normals)
	    hash.appendVec3(normal);
	hash.appendU32(
	    static_cast<uint32_t>(mesh.indices.size()));
	for (uint32_t index : mesh.indices)
	    hash.appendU32(index);
	hash.appendBox(mesh.bounds);
	const bool progressiveMesh = mesh.isProgressive();
	hash.appendByte(progressiveMesh ? 1 : 0);
	if (progressiveMesh) {
	    hash.appendU32(static_cast<uint32_t>(mesh.progressiveCuts.size()));
	    for (const Obol::ProgressiveTriangleCut &cut :
		 mesh.progressiveCuts) {
		hash.appendU32(cut.indexCount);
		hash.appendU32(cut.positionCount);
		hash.appendQuantization(cut.quantization);
	    }
	    hash.appendByte(mesh.progressiveMinimumCut);
	    hash.appendByte(mesh.progressiveResidentCut);
	    hash.appendVec3(mesh.progressiveQuantizationMinimum);
	    hash.appendVec3(mesh.progressiveQuantizationMaximum);
	    hash.appendU64(mesh.progressiveLineage);
	    hash.appendU32(static_cast<uint32_t>(
		mesh.progressiveClusters.size()));
	    hash.appendU32(mesh.progressiveClusterGridResolution);
	    for (const Obol::ProgressiveTriangleCluster &cluster :
		 mesh.progressiveClusters) {
		hash.appendBox(cluster.bounds);
		hash.appendByte(cluster.residentCut);
		hash.appendU32(static_cast<uint32_t>(cluster.ranges.size()));
		for (const Obol::ProgressiveTriangleClusterRange &range :
		     cluster.ranges) {
		    hash.appendU32(range.firstIndex);
		    hash.appendU32(range.indexCount);
		    hash.appendByte(range.activationCut);
		}
	    }
	}
    }
    const Obol::PartId contentId = hash.id();
    char digest[96] = {0};
    snprintf(digest, sizeof(digest), "%s:%016" PRIx64 "%016" PRIx64,
	kind, contentId.w1, contentId.w0);
    key.assign(digest);
    return 1;
}

static int
cad_vlist_part_geometry_supported(const SoBRLVListShape *shape,
				  const SoBRLVListShape **geomOut,
				  int *countOut)
{
    if (!shape || shape->editEmphasis.getValue() ||
	shape->selectedPrimitive.getNum() > 0 ||
	shape->highlightedPrimitive.getNum() > 0)
	return 0;

    const SoBRLVListShape *geom = shape->getGeometrySource();
    if (!geom)
	return 0;

    int n = geom->point.getNum();
    if (geom->command.getNum() < n)
	n = geom->command.getNum();
    if (n <= 0)
	return 0;

    if (geomOut)
	*geomOut = geom;
    if (countOut)
	*countOut = n;
    return 1;
}

static int
cad_vlist_part_geometry(const SoBRLVListShape *shape,
			Obol::PartGeometryBuilder &geometry)
{
    const SoBRLVListShape *geom = NULL;
    int n = 0;
    if (!cad_vlist_part_geometry_supported(shape, &geom, &n))
	return 0;

    Obol::WireRep wire;
    wire.bounds.makeEmpty();
    wire.segmentPoints.reserve(static_cast<size_t>(n) * 2u);
    wire.segmentIds.reserve(static_cast<size_t>(n));
    Obol::PointRep points;
    points.bounds.makeEmpty();
    points.positions.reserve(static_cast<size_t>(n));
    points.pointIds.reserve(static_cast<size_t>(n));
    SbBool haveLast = FALSE;
    int lastIndex = -1;
    uint32_t segmentIndex = 0;
    for (int i = 0; i < n; i++) {
	const int command = geom->command[i];
	if (command == SoBRLVListShape::POINT) {
	    SbVec3f point;
	    if (!cad_vlist_point(shape, geom, i, point))
		return 0;
	    points.positions.push_back(point);
	    points.pointIds.push_back(static_cast<uint32_t>(i));
	    points.bounds.extendBy(point);
	    SbColor color;
	    const SbBool colorValid = shape->getPointColor(i, color);
	    points.colorValid.push_back(colorValid ? 1u : 0u);
	    points.colors.push_back(colorValid ? color :
		SbColor(1.0f, 1.0f, 1.0f));
	    float scale = 0.0f;
	    const SbBool scaleValid = shape->getPointScale(i, scale);
	    points.scaleValid.push_back(scaleValid ? 1u : 0u);
	    points.scales.push_back(scaleValid ? scale : 0.0f);
	    if (scaleValid && scale > 0.0f) {
		const SbVec3f extent(scale, scale, scale);
		points.bounds.extendBy(point - extent);
		points.bounds.extendBy(point + extent);
	    }
	    SbVec3f normal;
	    const SbBool normalValid = shape->getPointNormal(i, normal);
	    points.normalValid.push_back(normalValid ? 1u : 0u);
	    points.normals.push_back(normalValid ? normal :
		SbVec3f(0.0f, 0.0f, 1.0f));
	    continue;
	}
	if (command == SoBRLVListShape::MOVE) {
	    haveLast = TRUE;
	    lastIndex = i;
	    continue;
	}
	if (command != SoBRLVListShape::DRAW || !haveLast ||
	    lastIndex < 0)
	    continue;

	SbVec3f a;
	SbVec3f b;
	if (!cad_vlist_point(shape, geom, lastIndex, a) ||
	    !cad_vlist_point(shape, geom, i, b))
	    return 0;
	wire.bounds.extendBy(a);
	wire.bounds.extendBy(b);
	wire.segmentPoints.push_back(a);
	wire.segmentPoints.push_back(b);
	wire.segmentIds.push_back(segmentIndex++);
	lastIndex = i;
    }

    if (wire.segmentPoints.empty() && points.positions.empty())
	return 0;
    if (!points.positions.empty())
	geometry.points = std::move(points);
    if (!wire.segmentPoints.empty())
	geometry.wire = std::move(wire);
    const char *source_type = shape->sourceType.getValue().getString();
    const char *geometry_kind = shape->geometryKind.getValue().getString();
    if (source_type && geometry_kind &&
	BU_STR_EQUAL(source_type, "proxy") &&
	(BU_STR_EQUAL(geometry_kind, "aabb") ||
	 BU_STR_EQUAL(geometry_kind, "obb"))) {
	/* The source remains a full conservative proxy for bounds and picks.
	 * SoCADAssembly alone decides whether its projected extent can be
	 * represented by one depth-tested pixel in this view. */
	geometry.subpixelProxyEligible = true;
	geometry.structuralProxy = true;
    }
    return 1;
}

static int
cad_mesh_part_geometry(const SoBRLMeshShape *shape,
		       Obol::PartGeometryBuilder &geometry)
{
    if (!shape || shape->editEmphasis.getValue() ||
	shape->selectedPrimitive.getNum() > 0 ||
	shape->highlightedPrimitive.getNum() > 0)
	return 0;

    const SoBRLMeshShape *geom = shape->getGeometrySource();
    if (!geom || geom->point.getNum() <= 0 ||
	geom->coordIndex.getNum() <= 0)
	return 0;

    Obol::TriMesh mesh;
    mesh.bounds.makeEmpty();
    mesh.positions.reserve(static_cast<size_t>(geom->point.getNum()));
    for (int i = 0; i < geom->point.getNum(); i++) {
	mesh.positions.push_back(geom->point[i]);
	mesh.bounds.extendBy(geom->point[i]);
    }
    mesh.indices.reserve(static_cast<size_t>(geom->coordIndex.getNum()));
    std::vector<int32_t> cornerIndices;
    cornerIndices.reserve(static_cast<size_t>(geom->coordIndex.getNum()));
    for (int i = 0; i < geom->coordIndex.getNum(); i++) {
	const int idx = geom->coordIndex[i];
	if (idx < 0 || idx >= geom->point.getNum())
	    return 0;
	mesh.indices.push_back(static_cast<uint32_t>(idx));
	cornerIndices.push_back(idx);
    }
    if (mesh.indices.empty() || mesh.bounds.isEmpty())
	return 0;

    std::vector<SbVec3f> cornerNormals;
    if (geom->normal.getNum() == geom->coordIndex.getNum()) {
	cornerNormals.reserve(static_cast<size_t>(geom->normal.getNum()));
	for (int i = 0; i < geom->normal.getNum(); ++i)
	    cornerNormals.push_back(geom->normal[i]);
    }
    sanitize_triangle_normals(cornerNormals, mesh.positions, cornerIndices);
    if (!canonicalize_corner_normal_mesh(mesh, cornerNormals))
	return 0;

    geometry.shaded = std::move(mesh);
    if (shape->hiddenLine.getValue() ||
	shape->drawMode.getValue() == BOBOL_LOD_DRAW_HIDDEN_LINE)
	(void)cad_mesh_append_hidden_line_edges(geometry);
    /* Legacy/evaluated CAD mesh records do not necessarily have PoP data,
     * but their complete vertex bounds are still conservative.  Authorize
     * the retained assembly to replace the whole occurrence with one
     * depth-tested point when its projected extent falls below the active
     * screen-error threshold.  This is the non-progressive escape path for
     * software rendering under frame pressure; selected occurrences are
     * promoted by SoCADAssembly and never remain collapsed. */
    geometry.subpixelProxyEligible = true;
    return 1;
}

std::string
cad_instance_key(const SoBRLDatabaseSource *source,
	const char *path, int ordinal)
{
    std::string key = source ?
	source_effective_instance_key(source).getString() : "";
    key.append("|");
    key.append(path ? path : "");
    key.append("#");
    char buf[64] = {0};
    snprintf(buf, sizeof(buf), "%d", ordinal);
    key.append(buf);
    return key;
}

SbMatrix
cad_instance_matrix(const SoBRLDatabaseSource *source,
		    const SbMatrix &localMatrix)
{
    SbMatrix matrix = localMatrix;
    if (source && source->drawMatrixValid.getValue())
	matrix.multRight(source->drawMatrix.getValue());
    canonicalize_affine_tail(matrix);
    return matrix;
}

SbMatrix
compact_mesh_asset_matrix(const SoBRLDatabaseSource *source,
	const BObolCompactInstanceEntry &entry)
{
    SbMatrix matrix = entry.sourceMeshRequest.meshAssetTransform;
    matrix.multRight(entry.placementTransform);
    return cad_instance_matrix(source, matrix);
}

struct cad_pending_part {
    Obol::PartId part;
    Obol::PartGeometryBuilder geometry;
};

struct cad_build_data {
    SoBRLDatabaseSource *source;
    std::map<std::string, Obol::PartId> partIdByKey;
    std::vector<cad_pending_part> parts;
    std::vector<Obol::InstanceUpdate> instances;
    std::vector<std::pair<Obol::InstanceId,
	SoBRLCadAssembly::InstanceSemantic>> semantics;
    std::vector<Obol::InstanceId> hiddenInstances;
    std::vector<Obol::InstanceId> selectedInstances;
    std::vector<Obol::InstanceId> unpickableInstances;
    int ordinal;
    int unsupported;
    int wireCount;
    int shadedCount;
};

static void
cad_add_part_if_needed(cad_build_data &data,
		       const std::string &partKey,
		       const Obol::PartGeometryBuilder &geometry,
		       Obol::PartId &partId)
{
    std::map<std::string, Obol::PartId>::iterator found =
	data.partIdByKey.find(partKey);
    if (found != data.partIdByKey.end()) {
	partId = found->second;
	return;
    }

    partId = Obol::CadIdBuilder::partId(partKey);
    data.partIdByKey[partKey] = partId;
    cad_pending_part update;
    update.part = partId;
    update.geometry = geometry;
    data.parts.push_back(update);
}

static SoBRLCadAssembly::InstanceSemantic
cad_vlist_semantic(const SoBRLVListShape *shape)
{
    SoBRLCadAssembly::InstanceSemantic semantic;
    if (!shape)
	return semantic;

    semantic.path = shape->sourcePath.getValue();
    semantic.sourceInstanceKey = shape->ownerSourceInstanceKey.getValue();
    semantic.sourceName = shape->sourceName.getValue();
    semantic.sourceType = shape->sourceType.getValue();
    semantic.sourceId = shape->sourceId.getValue();
    semantic.regionId = shape->regionId.getValue();
    semantic.airCode = shape->airCode.getValue();
    semantic.materialId = shape->materialId.getValue();
    semantic.los = shape->los.getValue();
    semantic.materialColorValid = shape->materialColorValid.getValue();
    semantic.materialColor = shape->materialColor.getValue();
    semantic.materialShader = shape->materialShader.getValue();
    semantic.editIntentId = shape->editIntentId.getValue();
    semantic.editIntentRole = shape->editIntentRole.getValue();
    semantic.primitiveKind = SoBRLPickDetail::LINE_SEGMENT;
    return semantic;
}

static SoBRLCadAssembly::InstanceSemantic
cad_mesh_semantic(const SoBRLMeshShape *shape)
{
    SoBRLCadAssembly::InstanceSemantic semantic;
    if (!shape)
	return semantic;

    semantic.path = shape->sourcePath.getValue();
    semantic.sourceInstanceKey = shape->ownerSourceInstanceKey.getValue();
    semantic.sourceName = shape->sourceName.getValue();
    semantic.sourceType = shape->sourceType.getValue();
    semantic.sourceId = shape->sourceId.getValue();
    semantic.regionId = shape->regionId.getValue();
    semantic.airCode = shape->airCode.getValue();
    semantic.materialId = shape->materialId.getValue();
    semantic.los = shape->los.getValue();
    semantic.materialColorValid = shape->materialColorValid.getValue();
    semantic.materialColor = shape->materialColor.getValue();
    semantic.materialShader = shape->materialShader.getValue();
    semantic.editIntentId = shape->editIntentId.getValue();
    semantic.editIntentRole = shape->editIntentRole.getValue();
    semantic.primitiveKind = SoBRLPickDetail::FACE;
    return semantic;
}

Obol::InstanceStyle
cad_source_style(const SoBRLDatabaseSource *source)
{
    Obol::InstanceStyle style;
    if (!source)
	return style;

    const SbColor materialColor =
	source->materialColorValid.getValue() ?
	source->materialColor.getValue() :
	(source->databaseMaterialColorValid.getValue() ?
	 source->databaseMaterialColor.getValue() :
	 source->color.getValue());

    style.hasColorOverride = true;
    cad_shape_color(source->selected.getValue(),
		    source->selectedColor.getValue(),
		    source->highlighted.getValue(),
		    source->highlightedColor.getValue(),
		    FALSE, source->ghostedColor.getValue(),
		    source->colorOverride.getValue(), source->color.getValue(),
		    TRUE, materialColor, source->color.getValue(),
		    source->transparency.getValue(), style.color);
    style.lineWidth = source->lineWidth.getValue() > 0 ?
		      static_cast<float>(source->lineWidth.getValue()) : 1.0f;
    if (source->lineStyle.getValue() != 0)
	style.linePattern = 0xcf33u;
    return style;
}

const char *
cad_source_leaf_name(const SoBRLDatabaseSource *source)
{
    if (!source)
	return "";
    const char *path = source->path.getValue().getString();
    if (!path || !path[0])
	return "";
    const char *slash = strrchr(path, '/');
    return (slash && slash[1]) ? slash + 1 : path;
}

Obol::InstanceId
cad_source_parent_instance(const SoBRLDatabaseSource *source)
{
    if (!source)
	return Obol::CadIdBuilder::rootInstance();
    const char *key = source->parentInstanceKey.getValue().getString();
    if (!key || !key[0])
	return Obol::CadIdBuilder::rootInstance();
    return Obol::CadIdBuilder::instanceId(key);
}

uint8_t
cad_source_boolean_operation(const SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;
    const int operation = source->booleanOperation.getValue();
    if (operation == SoBRLDatabaseSource::BOOLEAN_SUBTRACT)
	return 1;
    if (operation == SoBRLDatabaseSource::BOOLEAN_INTERSECT)
	return 2;
    return 0;
}

SoBRLCadAssembly::InstanceSemantic
cad_source_semantic(const SoBRLDatabaseSource *source,
		    SoBRLPickDetail::PrimitiveKind primitiveKind)
{
    SoBRLCadAssembly::InstanceSemantic semantic;
    if (!source)
	return semantic;

    semantic.path = source->path.getValue();
    semantic.sourceInstanceKey = source->instanceKey.getValue();
    semantic.sourceName = cad_source_leaf_name(source);
    semantic.sourceType =
	(primitiveKind == SoBRLPickDetail::FACE) ? "mesh-lod" : "proxy-lod";
    semantic.sourceId = source->sourceRevision.getValue();
    semantic.regionId = source->databaseRegionId.getValue();
    semantic.airCode = source->databaseAirCode.getValue();
    semantic.materialId = source->databaseMaterialId.getValue();
    semantic.los = source->databaseLos.getValue();
    semantic.materialColorValid =
	source->materialColorValid.getValue() ||
	source->databaseMaterialColorValid.getValue();
    semantic.materialColor =
	source->materialColorValid.getValue() ?
	source->materialColor.getValue() : source->databaseMaterialColor.getValue();
    semantic.materialShader = source->databaseMaterialShader.getValue();
    semantic.primitiveKind = primitiveKind;
    return semantic;
}

static void
cad_add_vlist_instance(cad_build_data &data,
		       SoBRLVListShape *shape,
		       const SbMatrix &localMatrix)
{
    if (!shape)
	return;
    Obol::PartGeometryBuilder geometry;
    if (!cad_vlist_part_geometry(shape, geometry)) {
	data.unsupported = 1;
	return;
    }

    std::string partKey;
    if (!cad_part_key_for_geometry("wire", geometry, partKey)) {
	data.unsupported = 1;
	return;
    }
    Obol::PartId partId;
    cad_add_part_if_needed(data, partKey, geometry, partId);

    const SbMatrix matrix = cad_instance_matrix(data.source, localMatrix);
    const std::string instanceKey =
	cad_instance_key(data.source, shape->sourcePath.getValue().getString(),
			 data.ordinal++);
    Obol::InstanceId instanceId =
	Obol::CadIdBuilder::instanceId(instanceKey);
    Obol::InstanceRecord record;
    record.part = partId;
    record.localToRoot = matrix;
    record.parent = cad_source_parent_instance(data.source);
    record.childName = shape->sourcePath.getValue().getString();
    record.occurrenceIndex = data.source->occurrenceIndex.getValue();
    record.boolOp = cad_source_boolean_operation(data.source);
    record.style = cad_vlist_style(shape);

    Obol::InstanceUpdate update;
    update.instance = instanceId;
    update.record = record;
    data.instances.push_back(update);
    if (!shape->visible.getValue())
	data.hiddenInstances.push_back(instanceId);
    if (shape->selected.getValue())
	data.selectedInstances.push_back(instanceId);
    if (!shape->selectable.getValue())
	data.unpickableInstances.push_back(instanceId);
    data.semantics.emplace_back(instanceId, cad_vlist_semantic(shape));
    data.wireCount++;
}

static void
cad_add_mesh_instance(cad_build_data &data,
		      SoBRLMeshShape *shape,
		      const SbMatrix &localMatrix)
{
    if (!shape)
	return;
    Obol::PartGeometryBuilder geometry;
    if (!cad_mesh_part_geometry(shape, geometry)) {
	data.unsupported = 1;
	return;
    }

    std::string partKey;
    if (!cad_part_key_for_geometry("mesh", geometry, partKey)) {
	data.unsupported = 1;
	return;
    }
    Obol::PartId partId;
    cad_add_part_if_needed(data, partKey, geometry, partId);

    const SbMatrix matrix = cad_instance_matrix(data.source, localMatrix);
    const std::string instanceKey =
	cad_instance_key(data.source, shape->sourcePath.getValue().getString(),
			 data.ordinal++);
    Obol::InstanceId instanceId =
	Obol::CadIdBuilder::instanceId(instanceKey);
    Obol::InstanceRecord record;
    record.part = partId;
    record.localToRoot = matrix;
    record.parent = cad_source_parent_instance(data.source);
    record.childName = shape->sourcePath.getValue().getString();
    record.occurrenceIndex = data.source->occurrenceIndex.getValue();
    record.boolOp = cad_source_boolean_operation(data.source);
    record.style = cad_mesh_style(shape);

    Obol::InstanceUpdate update;
    update.instance = instanceId;
    update.record = record;
    data.instances.push_back(update);
    if (!shape->visible.getValue())
	data.hiddenInstances.push_back(instanceId);
    if (shape->selected.getValue())
	data.selectedInstances.push_back(instanceId);
    if (!shape->selectable.getValue())
	data.unpickableInstances.push_back(instanceId);
    data.semantics.emplace_back(instanceId, cad_mesh_semantic(shape));
    data.shadedCount++;
}

static void
cad_collect_realized_node(cad_build_data &data,
			  SoNode *node,
			  const SbMatrix &matrix)
{
    if (!node || data.unsupported)
	return;

    if (node_is_source_placement_transform(node) ||
	node_is_auxiliary_vlist(node) || node_is_auxiliary_source(node))
	return;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	cad_add_vlist_instance(data, static_cast<SoBRLVListShape *>(node),
			       matrix);
	return;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	cad_add_mesh_instance(data, static_cast<SoBRLMeshShape *>(node),
			      matrix);
	return;
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;

    SoGroup *group = static_cast<SoGroup *>(node);
    SbMatrix childMatrix = matrix;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *child = group->getChild(i);
	if (!child)
	    continue;
	if (child->isOfType(SoMatrixTransform::getClassTypeId())) {
	    SoMatrixTransform *transform =
		static_cast<SoMatrixTransform *>(child);
	    childMatrix.multRight(transform->matrix.getValue());
	    continue;
	}
	cad_collect_realized_node(data, child, childMatrix);
    }
}

/* Only newly inserted, still-unpublished slots belong to this guard. Existing
 * values are read during preparation and written only in the no-allocation
 * commit. Capacity growth may survive a failed preparation; logical state may
 * not. Each occurrence uses a fixed number of these stack guards. */
template <typename Map>
class CompactPublicationSlot {
public:
    CompactPublicationSlot(Map &map, const typename Map::key_type &key,
	bool enabled = true) : map_(map), position_(map.end())
    {
	if (!enabled)
	    return;
	const auto result = map_.try_emplace(key);
	position_ = result.first;
	inserted_ = result.second;
    }
    ~CompactPublicationSlot()
    {
	if (inserted_)
	    map_.erase(position_);
    }
    CompactPublicationSlot(const CompactPublicationSlot &) = delete;
    CompactPublicationSlot &operator=(const CompactPublicationSlot &) = delete;
    typename Map::mapped_type &value() { return position_->second; }
    void commit() noexcept { inserted_ = false; }

private:
    Map &map_;
    typename Map::iterator position_;
    bool inserted_ = false;
};

template <typename T>
static void
compact_prepare_push(std::vector<T> &values)
{
    if (values.size() < values.capacity())
	return;
    /* Preserve amortized growth when the producer has not reserved a certified
     * population. Reserving size()+1 for every leaf makes streaming quadratic. */
    const size_t growthFactor = 2;
    const size_t capacity = values.size() <= values.max_size() / growthFactor ?
	std::max(size_t(1), values.size() * growthFactor) : values.max_size();
    if (values.size() == values.max_size())
	throw std::length_error("compact occurrence capacity exhausted");
    values.reserve(capacity);
}

template <typename Geometry>
static const char *
compact_geometry_kind(const Geometry &geometry)
{
    return geometry.shaded ? "mesh" : (geometry.points && !geometry.wire ? "point" : "wire");
}

class CompactPartPublication {
public:
    CompactPartPublication(BObolCompactInstanceIndex &index,
	const std::shared_ptr<const Obol::PartGeometry> &geometry,
	const char *overviewPath = nullptr,
	const BObolCompactInstanceEntry *previousOverview = nullptr) :
	index_(index), geometry_(geometry)
    {
	if (previousOverview) {
	    part_ = previousOverview->part;
	    for (BObolCompactPartReference &part : index_.parts) {
		if (part.part == part_) {
		    overviewPart_ = &part;
		    break;
		}
	    }
	    if (!overviewPart_)
		return;
	    previousGeometry_ = previousOverview->geometry.get();
	} else {
	    const auto cached = index_.partIdByGeometry.find(geometry.get());
	    if (!overviewPath && cached != index_.partIdByGeometry.end() &&
		bobol_compact_geometry_identity_matches(cached->second, geometry)) {
		part_ = cached->second.part;
		valid_ = true;
		return;
	    }
	    std::string partKey;
	    if (overviewPath) {
		/* Evolving overview extents never share the leaf proxy's part. */
		partKey = "compact-lod-overview:";
		partKey += overviewPath;
	    } else {
		const char *kind = compact_geometry_kind(*geometry);
		if (!cad_part_key_for_geometry(kind, *geometry, partKey))
		    return;
	    }
	    const auto existing = index_.partIdByKey.find(partKey);
	    appendPart_ = existing == index_.partIdByKey.end();
	    part_ = appendPart_ ? Obol::CadIdBuilder::partId(partKey) :
		existing->second;
	    if (appendPart_) {
		partSlot_.emplace(index_.partIdByKey, partKey);
		compact_prepare_push(index_.parts);
	    }
	}
	geometrySlot_.emplace(index_.partIdByGeometry, geometry.get());
	valid_ = true;
    }

    bool valid() const { return valid_; }
    const Obol::PartId &part() const { return part_; }
    void commit() noexcept
    {
	if (appendPart_) {
	    index_.parts.push_back({part_, geometry_});
	    partSlot_->value() = part_;
	    partSlot_->commit();
	}
	if (overviewPart_)
	    overviewPart_->geometry = geometry_;
	if (geometrySlot_) {
	    geometrySlot_->value() = {geometry_, part_};
	    geometrySlot_->commit();
	}
	if (previousGeometry_ && previousGeometry_ != geometry_.get()) {
	    const auto old = index_.partIdByGeometry.find(previousGeometry_);
	    if (old != index_.partIdByGeometry.end() && old->second.part == part_)
		index_.partIdByGeometry.erase(old);
	}
    }

private:
    BObolCompactInstanceIndex &index_;
    const std::shared_ptr<const Obol::PartGeometry> &geometry_;
    Obol::PartId part_;
    BObolCompactPartReference *overviewPart_ = nullptr;
    const Obol::PartGeometry *previousGeometry_ = nullptr;
    std::optional<CompactPublicationSlot<decltype(index_.partIdByKey)>> partSlot_;
    std::optional<CompactPublicationSlot<decltype(index_.partIdByGeometry)>>
	geometrySlot_;
    bool appendPart_ = false;
    bool valid_ = false;
};

static void
compact_index_count_entry(BObolCompactInstanceIndex &index,
	const BObolCompactInstanceEntry &entry, int contribution) noexcept
{
    index.shadedCount += entry.meshGeometry ? contribution : 0;
    index.wireCount += entry.wireGeometry || entry.pointGeometry ? contribution : 0;
    const bool progressive =
	bobol_compact_geometry_is_resident_progressive(entry.geometry);
    const auto count = [contribution](size_t &value, bool present) {
	if (present) {
	    if (contribution > 0)
		++value;
	    else
		--value;
	}
    };
    count(index.viewDependentCsgGeometryCount, entry.viewDependentCsgGeometry);
    count(index.sourceMeshRequestCount, entry.sourceMeshRequestValid);
    count(index.residentProgressiveGeometryCount, progressive);
    count(index.displayLodTargetCount, entry.sourceMeshRequestValid || progressive);
}

static void compact_index_bounds_add(BObolCompactInstanceIndex &index,
	BObolCompactInstanceEntry &entry);
static void compact_append_prepared(BObolCompactInstanceIndex &index,
	BObolCompactInstanceEntry &&entry, Obol::InstanceUpdate &&update,
	bool publishLookups);

static bool compact_sync_entry_from_source(BObolCompactInstanceEntry &entry,
	const SoBRLDatabaseSource *source);

static SoBRLCadAssembly::InstanceSemantic
compact_semantic_from_summary(const BObolRealizedShapeSummary &summary)
{
    SoBRLCadAssembly::InstanceSemantic semantic;
    semantic.path = summary.path;
    semantic.sourceName = summary.sourceName;
    semantic.sourceType = summary.sourceType;
    semantic.sourceId = summary.sourceId;
    semantic.regionId = summary.regionId;
    semantic.airCode = summary.airCode;
    semantic.materialId = summary.materialId;
    semantic.los = summary.los;
    semantic.materialColorValid = summary.materialColorValid;
    semantic.materialColor = summary.materialColor;
    semantic.materialShader = summary.materialShader;
    semantic.editIntentId = summary.editIntentId;
    semantic.editIntentRole = summary.editIntentRole;
    const char *geometryKind = summary.geometryKind.getString();
    if (geometryKind && strstr(geometryKind, "point"))
	semantic.primitiveKind = SoBRLPickDetail::POINT;
    else
	semantic.primitiveKind = summary.shapeKind ==
	    BObolRealizedShapeSummary::SHAPE_MESH ? SoBRLPickDetail::FACE :
	    SoBRLPickDetail::LINE_SEGMENT;
    return semantic;
}

static void
compact_add_occurrence(SoBRLDatabaseSource *source,
	BObolCompactInstanceIndex &index,
	const compact_occurrence_build &input,
	int &ordinal,
	int &unsupported,
	const BObolCompactOccurrenceRegistryState *presentation)
{
    const std::shared_ptr<const Obol::PartGeometry> &geometry =
	input.occurrence.geometry;
    if (!source || !geometry) {
	unsupported = 1;
	return;
    }

    const char *path = input.occurrence.summary.path.getString();
    if (!path || !path[0])
	path = source->path.getValue().getString();
    const bool overview = BU_STR_EQUAL(
	input.occurrence.summary.recordRole.getString(), "lod-overview");
    CompactPartPublication part(index, geometry, overview ? path : nullptr);
    if (!part.valid()) {
	unsupported = 1;
	return;
    }
    const Obol::PartId partId = part.part();

    SbMatrix geometryToSource = input.occurrence.geometryTransform;
    geometryToSource.multRight(input.occurrence.localTransform);
    const SbMatrix matrix = cad_instance_matrix(source, geometryToSource);
    const std::string instanceKey = cad_instance_key(source, path, ordinal);
    const Obol::InstanceId instanceId =
	Obol::CadIdBuilder::instanceId(instanceKey);

    BObolCompactInstanceEntry entry;
    entry.instance = instanceId;
    entry.part = partId;
    entry.geometry = geometry;
    entry.wireGeometry = geometry->wire ? TRUE : FALSE;
    entry.pointGeometry = geometry->points ? TRUE : FALSE;
    /* These flags describe resident draw channels, not source semantics.  In
     * particular, a source-backed progressive mesh begins with a wire AABB
     * and must select a wire draw path until a shaded payload is available. */
    entry.meshGeometry = geometry->shaded ? TRUE : FALSE;
    entry.viewDependentCsgGeometry =
	input.occurrence.viewDependentCsgGeometry;
    entry.lodBacked = input.occurrence.lodBacked;
    entry.sourceMeshRequestValid = input.occurrence.sourceMeshRequestValid;
    if (entry.sourceMeshRequestValid) {
	entry.sourceMeshRequest = input.occurrence.sourceMeshRequest;
    }
    entry.localToSource = matrix;
    entry.geometryTransform = input.occurrence.geometryTransform;
    entry.placementTransform = input.occurrence.localTransform;
    entry.localTransform = geometryToSource;
    entry.semantic = input.semantic;
    if (entry.semantic.path.getLength() == 0)
	entry.semantic.path = path;
    entry.instanceKey = instanceKey.c_str();
    entry.semantic.sourceInstanceKey = instanceKey.c_str();
    entry.authoredVisible = input.occurrence.summary.visible;
    entry.visible = entry.authoredVisible;
    entry.selectable = input.occurrence.summary.selectable;
    entry.selected = input.occurrence.summary.selected;
    entry.authoredHighlighted = input.occurrence.summary.highlighted;
    entry.highlighted = compact_effective_highlight(entry);
    entry.shapeSummary = input.occurrence.summary;
    entry.occurrenceIndex = input.occurrence.occurrenceIndex;
    entry.booleanOperation = input.occurrence.booleanOperation;
    compact_sync_entry_from_source(entry, source);
    if (input.stylesValid) {
	entry.normalStyle = input.normalStyle;
	entry.selectedStyle = input.selectedStyle;
	entry.highlightedStyle = input.highlightedStyle;
    }
    if (input.dashed ||
	entry.booleanOperation == SoBRLDatabaseSource::BOOLEAN_SUBTRACT) {
	constexpr uint16_t subtractLinePattern = 0xcf33u;
	entry.normalStyle.linePattern = subtractLinePattern;
	entry.selectedStyle.linePattern = subtractLinePattern;
	entry.highlightedStyle.linePattern = subtractLinePattern;
    }
    entry.style = compact_effective_style(entry);
    if (presentation)
	compact_prepare_occurrence_presentation(entry, *presentation);
    compact_sync_shape_summary(entry);

    Obol::InstanceRecord record;
    record.part = partId;
    record.localToRoot = matrix;
    record.parent = cad_source_parent_instance(source);
    record.childName = path;
    record.occurrenceIndex = entry.occurrenceIndex;
    record.boolOp = entry.booleanOperation ==
	SoBRLDatabaseSource::BOOLEAN_SUBTRACT ? 1 :
	(entry.booleanOperation == SoBRLDatabaseSource::BOOLEAN_INTERSECT ? 2 : 0);
    /* PartGeometry::structuralProxy keeps an overview wire-visible in shaded
     * mode.  InstanceRecord::lodStructuralProxy has the narrower meaning
     * "unresolved LoD leaf" and must not make the aggregate overview enter
     * leaf-repair or convergence accounting. */
    record.lodStructuralProxy = geometry->structuralProxy && entry.lodBacked;
    record.style = entry.style;

    Obol::InstanceUpdate update;
    update.instance = instanceId;
    update.record = record;
    compact_append_prepared(index, std::move(entry), std::move(update),
	presentation != nullptr);
    part.commit();
    ++ordinal;
}

static void
compact_index_bounds_add(BObolCompactInstanceIndex &index,
	BObolCompactInstanceEntry &entry)
{
    entry.sourceBounds = database_source_transform_bounds(
	compact_part_geometry_bounds(entry.geometry), entry.localTransform);
    if (entry.sourceBounds.isEmpty())
	return;

    index.sourceBounds.extendBy(entry.sourceBounds);
}

static bool
compact_bounds_removal_requires_rebuild(const SbBox3f &aggregate, const SbBox3f &entry)
{
    if (entry.isEmpty() || aggregate.isEmpty())
	return false;
    for (int axis = 0; axis < 3; ++axis)
	if (entry.getMin()[axis] <= aggregate.getMin()[axis] ||
	    entry.getMax()[axis] >= aggregate.getMax()[axis])
	    return true;
    return false;
}

static void
compact_index_bounds_remove(BObolCompactInstanceIndex &index,
	BObolCompactInstanceEntry &entry)
{
    index.sourceBoundsDirty = index.sourceBoundsDirty ||
	compact_bounds_removal_requires_rebuild(index.sourceBounds, entry.sourceBounds);
    entry.sourceBounds.makeEmpty();
}

static SbBool
compact_index_source_bounds(BObolCompactInstanceIndex &index,
	SbBox3f &bounds)
{
    if (index.sourceBoundsDirty) {
	index.sourceBounds.makeEmpty();
	for (const BObolCompactInstanceEntry &entry : index.entries)
	    if (!entry.sourceBounds.isEmpty())
		index.sourceBounds.extendBy(entry.sourceBounds);
	index.sourceBoundsDirty = false;
    }

    bounds = index.sourceBounds;
    if (bounds.isEmpty())
	return FALSE;
    return TRUE;
}

static void
compact_rebuild_entry_index(BObolCompactInstanceIndex &index)
{
    index.entryIndex.clear();
    index.entryIndex.reserve(index.entries.size());
    index.entryIndexByKey.clear();
    index.entryIndexByKey.reserve(index.entries.size());
    index.entryIndexByPath.clear();
    index.entryIndexByPath.reserve(index.entries.size());
    index.entryIndexByOrderedPath.clear();

    index.entryIndicesByLeaf.clear();
    index.entryIndicesByLeaf.reserve(index.entries.size());
    index.entryIndicesBySourceName.clear();
    index.entryIndicesBySourceName.reserve(index.entries.size());
    index.partReferenceCounts.clear();
    index.partReferenceCounts.reserve(index.parts.size());
    index.sourceBounds.makeEmpty();
    index.sourceBoundsDirty = false;
    for (size_t i = 0; i < index.entries.size(); i++) {
	index.entryIndex[index.entries[i].instance] = i;
	const SbString occurrenceKey =
	    compact_instance_identity(index.entries[i]);
	if (occurrenceKey.getLength() > 0)
	    index.entryIndexByKey[occurrenceKey.getString()] = i;
	const char *path = database_source_skip_leading_slash(
	    index.entries[i].semantic.path.getString());
	if (path && path[0]) {
	    index.entryIndexByPath[path] = i;
	    index.entryIndexByOrderedPath[path] = i;
	}
	const std::string leaf = database_source_leaf_component(
	    index.entries[i].semantic.path);
	if (!leaf.empty())
	    index.entryIndicesByLeaf[leaf].push_back(i);
	const char *sourceName = index.entries[i].semantic.sourceName.getString();
	if (sourceName && sourceName[0])
	    index.entryIndicesBySourceName[sourceName].push_back(i);
	index.partReferenceCounts[index.entries[i].part]++;
	compact_index_bounds_add(index, index.entries[i]);
    }
}

static bool
compact_replaces_overview_baseline(
    const BObolRealizedShapeSummary &previous,
    const BObolRealizedShapeSummary &replacement)
{
    return BU_STR_EQUAL(previous.recordRole.getString(), "lod-overview") &&
	!BU_STR_EQUAL(replacement.recordRole.getString(), "lod-overview");
}

static void
compact_apply_occurrence_baseline(BObolCompactInstanceEntry &entry,
    const BObolRealizedShapeSummary &summary)
{
    entry.authoredVisible = summary.visible;
    entry.visible = compact_effective_authored_visibility(entry);
    entry.selectable = summary.selectable;
    entry.selected = summary.selected;
    entry.authoredHighlighted = summary.highlighted;
    entry.highlighted = compact_effective_highlight(entry);
}

static void
compact_prepare_installation(const BObolCompactInstanceIndex *current,
    BObolCompactInstanceIndex &next,
    const BObolCompactOccurrenceRegistryState &presentation)
{
    std::unordered_map<Obol::InstanceId,
	const BObolCompactInstanceEntry *, std::hash<Obol::InstanceId>> old;
    std::unordered_map<std::string, const BObolCompactInstanceEntry *>
	oldByPath;
    if (current) {
	old.reserve(current->entries.size());
	for (const BObolCompactInstanceEntry &entry : current->entries) {
	    auto oldInserted = old.emplace(entry.instance, &entry);
	    if (!oldInserted.second)
		oldInserted.first->second = NULL;
	    std::string path = database_source_skip_leading_slash(
		entry.semantic.path.getString());
	    if (!path.empty()) {
		auto inserted = oldByPath.emplace(path, &entry);
		if (!inserted.second)
		    inserted.first->second = NULL;
	    }
	}
    }

    std::unordered_set<Obol::InstanceId, std::hash<Obol::InstanceId>>
	assignedInstances;
    assignedInstances.reserve(next.entries.size());
    for (size_t i = 0; i < next.entries.size(); i++) {
	BObolCompactInstanceEntry &entry = next.entries[i];
	auto found = old.find(entry.instance);
	const BObolCompactInstanceEntry *previousEntry =
	    found != old.end() ? found->second : NULL;
	if (!previousEntry) {
	    std::string path = database_source_skip_leading_slash(
		entry.semantic.path.getString());
	    const auto pathIt = oldByPath.find(path);
	    if (pathIt != oldByPath.end())
		previousEntry = pathIt->second;
	}
	if (previousEntry) {
	    const BObolCompactInstanceEntry &previous = *previousEntry;
	    const bool replacesOverview =
		compact_replaces_overview_baseline(previous.shapeSummary,
		    entry.shapeSummary);
	    /* A source name or list position is not an occurrence identity:
	     * multiple CAD paths routinely share both geometry and leaf names,
	     * and a progressive current index may contain only one of them.
	     * Transfer a retained handle only for an exact instance/path match,
	     * and never assign one handle to two authoritative entries. */
	    if (entry.instance != previous.instance &&
		assignedInstances.find(previous.instance) ==
		    assignedInstances.end()) {
		entry.instance = previous.instance;
		entry.instanceKey = previous.instanceKey;
		entry.semantic.sourceInstanceKey = previous.instanceKey;
		if (i < next.instances.size()) {
		    next.instances[i].instance = previous.instance;
		}
	    }
	    if (!replacesOverview) {
		entry.authoredVisible = previous.authoredVisible;
		entry.selectable = previous.selectable;
		entry.selected = previous.selected;
		entry.authoredHighlighted = previous.authoredHighlighted;
	    }
	    entry.presentationVisibleValid = previous.presentationVisibleValid;
	    entry.presentationVisible = previous.presentationVisible;
	    entry.presentationHighlightedValid =
		previous.presentationHighlightedValid;
	    entry.presentationHighlighted = previous.presentationHighlighted;
	    entry.presentationTransparencyValid =
		previous.presentationTransparencyValid;
	    entry.presentationTransparency = previous.presentationTransparency;
	    if (replacesOverview)
		compact_apply_occurrence_baseline(entry, entry.shapeSummary);
	    else {
		entry.visible = previous.visible;
		entry.highlighted = previous.highlighted;
	    }
	    entry.geometryRevision = previous.geometryRevision;
	    if (entry.part != previous.part ||
		entry.lodBacked != previous.lodBacked)
		entry.geometryRevision = compact_next_revision(
		    entry.geometryRevision);
	    entry.placementRevision = previous.placementRevision;
	    if (!entry.localToSource.equals(previous.localToSource,
		    0.000001f))
		entry.placementRevision = compact_next_revision(
		    entry.placementRevision);
	    entry.visibilityRevision = previous.visibilityRevision;
	    entry.selectionRevision = previous.selectionRevision;
	    entry.appearanceRevision = previous.appearanceRevision;
	    entry.semanticRevision = previous.semanticRevision;
	    if (!compact_semantic_equal(entry.semantic, previous.semantic))
		compact_note_semantic_change(entry);
	    if (replacesOverview &&
		(entry.visible != previous.visible ||
		 entry.selectable != previous.selectable))
		entry.visibilityRevision = compact_next_revision(
		    entry.visibilityRevision);
	    if (replacesOverview &&
		(entry.selected != previous.selected ||
		 entry.highlighted != previous.highlighted))
		entry.selectionRevision = compact_next_revision(
		    entry.selectionRevision);
	}
	assignedInstances.insert(entry.instance);
	entry.style = compact_effective_style(entry);
	if (previousEntry && !compact_appearance_equal(entry, *previousEntry))
	    entry.appearanceRevision = compact_next_revision(
		entry.appearanceRevision);
	if (i < next.instances.size())
	{
	    next.instances[i].record.style = entry.style;
	    next.instances[i].record.lodStructuralProxy =
		entry.geometry && entry.geometry->structuralProxy &&
		entry.lodBacked;
	}
	compact_sync_shape_summary(entry);
    }
    compact_rebuild_entry_index(next);
    compact_prepare_registry_presentation(next, presentation);
}

int
SoBRLDatabaseSource::setCompactOccurrence(
    const BObolCompactOccurrence &occurrence)
{
    if (!occurrence.geometry)
	return 0;

    BObolPerformanceTimer timer(BOBOL_PERF_CAD_COMPACT_US);
    if (timer.active())
	bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_ATTEMPTS, 1);

    auto next = std::make_unique<BObolCompactInstanceIndex>();
    compact_occurrence_build input;
    input.occurrence = occurrence;
    input.semantic = compact_semantic_from_summary(occurrence.summary);
    int ordinal = 0;
    int unsupported = 0;
    compact_add_occurrence(this, *next, input, ordinal, unsupported);
    if (unsupported || next->entries.empty())
	return 0;

    this->installCompactInstanceIndex(std::move(next), FALSE,
	this->d->compactHandleSourceId);
    remove_non_auxiliary_children(this);
    this->markCompiledAssemblyDirty();

    bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_SOURCES, 1);
    bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_INSTANCES,
	static_cast<uint64_t>(this->d->compactIndex->entries.size()));
    return static_cast<int>(this->d->compactIndex->entries.size());
}

int
SoBRLDatabaseSource::setCompactOccurrenceRegistry(
    const std::vector<BObolCompactOccurrence> &occurrences)
{
    if (occurrences.empty())
	return 0;

    BObolPerformanceTimer timer(BOBOL_PERF_CAD_COMPACT_US);
    if (timer.active())
	bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_ATTEMPTS, 1);

    auto next = std::make_unique<BObolCompactInstanceIndex>();
    int ordinal = 0;
    int unsupported = 0;
    for (const BObolCompactOccurrence &occurrence : occurrences) {
	if (!occurrence.geometry) {
	    unsupported = 1;
	    break;
	}
	compact_occurrence_build input;
	input.occurrence = occurrence;
	input.semantic = compact_semantic_from_summary(occurrence.summary);
	compact_add_occurrence(this, *next, input, ordinal, unsupported);
	if (unsupported)
	    break;
    }
    if (unsupported || next->entries.size() != occurrences.size())
	return 0;

    this->installCompactInstanceIndex(std::move(next), TRUE,
	this->d->compactHandleSourceId);
    remove_non_auxiliary_children(this);
    this->markCompiledAssemblyDirty();

    bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_SOURCES, 1);
    bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_INSTANCES,
	static_cast<uint64_t>(this->d->compactIndex->entries.size()));
    return this->d->compactIndex->entries.size() >
	static_cast<size_t>(INT_MAX) ? INT_MAX :
	static_cast<int>(this->d->compactIndex->entries.size());
}

/* The deque append is the last allocating operation. All parallel arrays,
 * lookup slots, memberships, counters and bounds then commit without failure.
 * Detached full-source construction defers lookups to its installation pass. */
static void
compact_append_prepared(BObolCompactInstanceIndex &index,
	BObolCompactInstanceEntry &&entry, Obol::InstanceUpdate &&update,
	bool publishLookups)
{
    static_assert(std::is_nothrow_move_constructible<
	BObolCompactInstanceEntry>::value, "compact commit must not allocate");
    static_assert(std::is_nothrow_move_constructible<
	Obol::InstanceUpdate>::value, "instance commit must not allocate");
    compact_prepare_push(index.instances);
    if (!entry.visible)
	compact_prepare_push(index.hiddenInstances);
    if (entry.selected)
	compact_prepare_push(index.selectedInstances);
    if (!entry.selectable)
	compact_prepare_push(index.unpickableInstances);

    const size_t entryIndex = index.entries.size();
    const std::string key = publishLookups ?
	compact_instance_identity(entry).getString() : "";
    const char *path = publishLookups ? database_source_skip_leading_slash(
	entry.semantic.path.getString()) : "";
    const std::string leaf = publishLookups ?
	database_source_leaf_component(entry.semantic.path) : "";
    const char *name = publishLookups ? entry.semantic.sourceName.getString() : "";
    CompactPublicationSlot idSlot(index.entryIndex, entry.instance, publishLookups);
    CompactPublicationSlot keySlot(index.entryIndexByKey, key, !key.empty());
    CompactPublicationSlot pathSlot(index.entryIndexByPath, std::string(path), path[0]);
    CompactPublicationSlot orderedSlot(index.entryIndexByOrderedPath,
	std::string(path), path[0]);
    CompactPublicationSlot leafSlot(index.entryIndicesByLeaf, leaf, !leaf.empty());
    CompactPublicationSlot nameSlot(index.entryIndicesBySourceName,
	std::string(name), name[0]);
    CompactPublicationSlot partSlot(index.partReferenceCounts, entry.part,
	publishLookups);
    if (!leaf.empty())
	compact_prepare_push(leafSlot.value());
    if (name[0])
	compact_prepare_push(nameSlot.value());

    const bool hasPath = path[0];
    const bool hasName = name[0];
    index.entries.push_back(std::move(entry));
    const BObolCompactInstanceEntry &published = index.entries.back();
    index.instances.push_back(std::move(update));
    if (!published.visible)
	index.hiddenInstances.push_back(published.instance);
    if (published.selected)
	index.selectedInstances.push_back(published.instance);
    if (!published.selectable)
	index.unpickableInstances.push_back(published.instance);
    if (publishLookups) {
	idSlot.value() = entryIndex;
	idSlot.commit();
	++partSlot.value();
	partSlot.commit();
	compact_index_bounds_add(index, index.entries.back());
    }
    if (!key.empty()) {
	keySlot.value() = entryIndex;
	keySlot.commit();
    }
    if (hasPath) {
	pathSlot.value() = entryIndex;
	pathSlot.commit();
	orderedSlot.value() = entryIndex;
	orderedSlot.commit();
    }
    if (!leaf.empty()) {
	leafSlot.value().push_back(entryIndex);
	leafSlot.commit();
    }
    if (hasName) {
	nameSlot.value().push_back(entryIndex);
	nameSlot.commit();
    }
    compact_index_count_entry(index, published, 1);
}

/* Drawing-data tier of a compact occurrence: a coarse proxy box is superseded by
 * a tighter proxy which is superseded by realized geometry.  Progressive
 * streaming upgrades a leaf's occurrence in place when a higher tier arrives. */
static int
compact_geometry_tier(const char *geometryKind)
{
    if (!geometryKind || !geometryKind[0])
	return 2;
    if (BU_STR_EQUAL(geometryKind, "aabb") ||
	BU_STR_EQUAL(geometryKind, "overview-aabb"))
	return 0;
    if (BU_STR_EQUAL(geometryKind, "obb"))
	return 1;
    return 2;
}

static void
compact_set_instance_membership(std::vector<Obol::InstanceId> &instances,
	const Obol::InstanceId &instance, bool member)
{
    const std::vector<Obol::InstanceId>::iterator found =
	std::find(instances.begin(), instances.end(), instance);
    if (member) {
	if (found == instances.end())
	    instances.push_back(instance);
	return;
    }
    if (found != instances.end())
	instances.erase(found);
}

/* The caller prepares the entry, part and reference-count slot first. Keep
 * the stable deque address while publishing all derived channels together. */
static void
compact_commit_entry_geometry(BObolCompactInstanceIndex &index, size_t position,
    BObolCompactInstanceEntry &&entry) noexcept
{
    static_assert(std::is_nothrow_move_constructible<BObolCompactInstanceEntry>::value,
	"compact replacement must not allocate");
    BObolCompactInstanceEntry &previous = index.entries[position];
    compact_index_bounds_remove(index, previous);
    --index.partReferenceCounts.find(previous.part)->second;
    compact_index_count_entry(index, previous, -1);
    previous.~BObolCompactInstanceEntry();
    new (&previous) BObolCompactInstanceEntry(std::move(entry));
    ++index.partReferenceCounts.find(previous.part)->second;
    compact_index_count_entry(index, previous, 1);
    compact_index_bounds_add(index, previous);
    Obol::InstanceUpdate &update = index.instances[position];
    update.record.part = previous.part;
    update.record.localToRoot = previous.localToSource;
    update.record.lodStructuralProxy = previous.geometry->structuralProxy && previous.lodBacked;
    update.record.style = previous.style;
}

/* Replace one leaf occurrence's drawing data in place: swap its part geometry,
 * transform, flags, and summary for the incoming higher-tier occurrence while
 * preserving the entry's instance identity, semantic key, and leaf runtime
 * state.  An overview-to-leaf transition adopts the leaf's authored baseline.
 * instances[i] runs parallel to entries[i] (both are appended together by
 * compact_add_occurrence and never reordered), so the instance record is
 * updated at the same index. */
static bool
compact_index_replace_entry_geometry(SoBRLDatabaseSource *source,
	BObolCompactInstanceIndex &index, size_t entryIdx,
	const BObolCompactOccurrence &occurrence,
	const BObolCompactOccurrenceRegistryState &presentation)
{
    const std::shared_ptr<const Obol::PartGeometry> &geometry =
	occurrence.geometry;
    if (!source || !geometry || entryIdx >= index.entries.size() ||
	entryIdx >= index.instances.size())
	return false;

    BObolCompactInstanceEntry &previous = index.entries[entryIdx];
    BObolCompactInstanceEntry entry = previous;
    const bool evolvingOverview =
	BU_STR_EQUAL(entry.shapeSummary.recordRole.getString(),
	    "lod-overview") &&
	BU_STR_EQUAL(occurrence.summary.recordRole.getString(),
	    "lod-overview");
    const bool replacesOverview = compact_replaces_overview_baseline(
	entry.shapeSummary, occurrence.summary);

    CompactPartPublication part(index, geometry, nullptr,
	evolvingOverview ? &previous : nullptr);
    if (!part.valid())
	return false;
    const Obol::PartId newPartId = part.part();

    SbMatrix geometryToSource = occurrence.geometryTransform;
    geometryToSource.multRight(occurrence.localTransform);
    const SbMatrix matrix = cad_instance_matrix(source, geometryToSource);

    entry.part = newPartId;
    entry.geometry = geometry;
    entry.wireGeometry = geometry->wire ? TRUE : FALSE;
    entry.pointGeometry = geometry->points ? TRUE : FALSE;
    entry.meshGeometry = geometry->shaded ? TRUE : FALSE;
    entry.viewDependentCsgGeometry = occurrence.viewDependentCsgGeometry;
    entry.lodBacked = occurrence.lodBacked;
    entry.sourceMeshRequestValid = occurrence.sourceMeshRequestValid;
    if (entry.sourceMeshRequestValid)
	entry.sourceMeshRequest = occurrence.sourceMeshRequest;
    else
	entry.sourceMeshRequest.clear();
    /* A stable overview part can change from normalized to source-local
     * coordinates. Its shared-part update does not publish instance placement;
     * advance that independent revision when the replacement changes it. */
    if (entry.localToSource != matrix)
	entry.placementRevision = compact_next_revision(entry.placementRevision);
    entry.localToSource = matrix;
    entry.geometryTransform = occurrence.geometryTransform;
    entry.placementTransform = occurrence.localTransform;
    entry.localTransform = geometryToSource;
    /* Keep the entry's path/sourceName/instanceKey (leaf identity) untouched;
     * only the drawing-data descriptors change. */
    entry.shapeSummary.shapeKind = occurrence.summary.shapeKind;
    entry.shapeSummary.geometryKind = occurrence.summary.geometryKind;
    entry.shapeSummary.recordRole = occurrence.summary.recordRole;
    entry.shapeSummary.sourceType = occurrence.summary.sourceType;
    if (replacesOverview) {
	/* A whole-target overview is intentionally not selectable.  When the
	 * draw root is itself a leaf, its authoritative occurrence has the same
	 * path and upgrades that overview in place.  Presentation overlays still
	 * belong to the path, but the replacement's authored interaction state is
	 * the new baseline; inheriting the overview baseline makes a fully drawn
	 * primitive visible yet impossible to select. */
	compact_apply_occurrence_baseline(entry, occurrence.summary);
	compact_prepare_occurrence_presentation(entry, presentation);
	entry.visibilityRevision = previous.visibilityRevision;
	entry.selectionRevision = previous.selectionRevision;
	if (entry.visible != previous.visible ||
	    entry.selectable != previous.selectable)
	    entry.visibilityRevision = compact_next_revision(
		entry.visibilityRevision);
	if (entry.selected != previous.selected ||
	    entry.highlighted != previous.highlighted)
	    entry.selectionRevision = compact_next_revision(
		entry.selectionRevision);
    }
    entry.style = compact_effective_style(entry);
    compact_sync_shape_summary(entry);
    /*
     * This revision is the sparse publication contract between the compact
     * occurrence registry and its retained SoCADAssembly presentation.
     * Without advancing it, the changed-entry journal is correctly delivered
     * but compactViewLodAssembly filters the entry as already presented:
     * metadata then reports the richer part while the renderer continues to
     * own the old structural box indefinitely.
     */
    entry.geometryRevision =
	compact_next_revision(entry.geometryRevision);

    CompactPublicationSlot newReference(index.partReferenceCounts, newPartId);
    if (!entry.visible && previous.visible)
	compact_prepare_push(index.hiddenInstances);
    if (!entry.selectable && previous.selectable)
	compact_prepare_push(index.unpickableInstances);
    if (entry.selected && !previous.selected)
	compact_prepare_push(index.selectedInstances);

    /* All preparation succeeded; the remaining publication cannot allocate. */
    part.commit();
    compact_commit_entry_geometry(index, entryIdx, std::move(entry));
    newReference.commit();
    if (replacesOverview) {
	compact_set_instance_membership(index.hiddenInstances, previous.instance,
	    !previous.visible);
	compact_set_instance_membership(index.unpickableInstances, previous.instance,
	    !previous.selectable);
	compact_set_instance_membership(index.selectedInstances, previous.instance,
	    previous.selected);
    }
    return true;
}

/* Enrich an already-published proxy with its source-backed PoP contract
 * without replacing the drawable proxy.  Cold coverage intentionally
 * publishes one shared unit box plus a per-occurrence normalization matrix.
 * The later authoritative tree walk may carry an independently allocated
 * exact AABB.  Replacing the shared box merely because the source request has
 * arrived defeats instancing, and (before proxy and asset transforms were
 * separated) also risked applying box normalization to the mesh asset.
 *
 * Geometry tier changes still go through
 * compact_index_replace_entry_geometry. */
static bool
compact_index_merge_source_contract(BObolCompactInstanceIndex &index,
	size_t entryIdx, const BObolCompactOccurrence &occurrence)
{
    if (entryIdx >= index.entries.size())
	return false;

    BObolCompactInstanceEntry &previous = index.entries[entryIdx];
    BObolCompactInstanceEntry entry = previous;
    bool changed = false;
    if (occurrence.viewDependentCsgGeometry &&
	!entry.viewDependentCsgGeometry) {
	entry.viewDependentCsgGeometry = TRUE;
	changed = true;
    }
    if (occurrence.lodBacked && !entry.lodBacked) {
	entry.lodBacked = TRUE;
	changed = true;
    }
    if (occurrence.sourceMeshRequestValid) {
	if (!entry.sourceMeshRequestValid) {
	    entry.sourceMeshRequestValid = TRUE;
	}
	entry.sourceMeshRequest = occurrence.sourceMeshRequest;
	changed = true;
    }
    if (!changed)
	return false;
    /* Source semantics become authoritative with the request, but the
     * currently drawn arrays, proxy transform, placement, part, and runtime
     * visibility/selection state remain untouched.  The unresolved-LoD role
     * is nevertheless an occurrence presentation property: it must reach the
     * retained assembly so convergence cannot mistake this box for authored
     * geometry.  Advance the existing sparse record revision and let Obol's
     * instance-attribute journal patch the one flag without rebuilding the
     * shared part or frame plan. */
    entry.geometryRevision =
	compact_next_revision(entry.geometryRevision);
    entry.shapeSummary.shapeKind = occurrence.summary.shapeKind;
    entry.shapeSummary.sourceType = occurrence.summary.sourceType;
    entry.shapeSummary.sourceId = occurrence.summary.sourceId;
    compact_sync_shape_summary(entry);
    compact_index_count_entry(index, previous, -1);
    previous.~BObolCompactInstanceEntry();
    new (&previous) BObolCompactInstanceEntry(std::move(entry));
    compact_index_count_entry(index, previous, 1);
    if (entryIdx < index.instances.size())
	index.instances[entryIdx].record.lodStructuralProxy =
	    previous.geometry && previous.geometry->structuralProxy &&
	    previous.lodBacked;
    return true;
}

static bool
compact_complete_stream_leaf_frontier(
    const BObolCompactInstanceIndex &index, size_t expectedLeafCount)
{
    /* The expected population excludes the temporary overview.  Avoid the
     * linear confirmation until the index is large enough to contain both
     * that overview and every expected leaf. */
    if (!expectedLeafCount || index.entries.size() <= expectedLeafCount)
	return false;

    size_t leafCount = 0;
    for (const BObolCompactInstanceEntry &entry : index.entries) {
	if (BU_STR_EQUAL(entry.shapeSummary.recordRole.getString(),
		"lod-overview"))
	    continue;
	if (++leafCount >= expectedLeafCount)
	    return true;
    }
    return false;
}

int
SoBRLDatabaseSource::mergeCompactOccurrences(
    const std::vector<BObolCompactOccurrence> &occurrences,
    SbBool authoritativeGeometry)
{
    return this->mergeCompactOccurrences(occurrences, authoritativeGeometry,
	nullptr, nullptr, nullptr, 0);
}

int
SoBRLDatabaseSource::mergeCompactOccurrences(
    const std::vector<BObolCompactOccurrence> &occurrences,
    SbBool authoritativeGeometry,
    PublicationCommit committed,
    void *context,
    SbBool *batchCompleted,
    size_t reserveCapacity)
{
    if (batchCompleted)
	*batchCompleted = FALSE;
    if (occurrences.empty())
	return 0;

    if (reserveCapacity &&
	this->reserveCompactStreamCapacity(reserveCapacity) < 0)
	return -1;
    if (!this->d->compactIndex) {
	this->d->compactIndex = new BObolCompactInstanceIndex;
	this->d->compactIndexActive = TRUE;
	this->d->compactOccurrenceRegistry = TRUE;
    }
    BObolCompactInstanceIndex &index = *this->d->compactIndex;

    const size_t firstNew = index.entries.size();
    int ordinal = static_cast<int>(index.entries.size());
    int changed = 0;
    bool mergedSourceContract = false;
    bool overviewArrived = false;
    std::vector<size_t> changedEntries;
    changedEntries.reserve(occurrences.size());
    const bool denyPartialMerge = occurrences.size() > 1 &&
	bobol_transaction_fault_requested(
	    BObolTransactionFaultPoint::SOURCE_STREAM_PARTIAL_MERGE_AFTER_COMPLETION);
    const auto publish = [&]() {
	if (!changed)
	    return;
	if (overviewArrived && this->d->compactOverviewState ==
		BObolCompactOccurrenceRegistryState::OverviewState::Absent)
	    this->d->compactOverviewState =
		BObolCompactOccurrenceRegistryState::OverviewState::Visible;
	if (this->d->compactOverviewState ==
		BObolCompactOccurrenceRegistryState::OverviewState::Visible &&
	    this->d->compactExpectedInstanceCountCertified &&
	    compact_complete_stream_leaf_frontier(index,
		this->d->compactExpectedInstanceCount))
	    this->d->compactOverviewState =
		BObolCompactOccurrenceRegistryState::OverviewState::RetirementPending;
	this->d->compactOccurrenceRegistry = TRUE;
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty(changedEntries);
	this->markDisplayMeshLodDirty(changedEntries,
	    index.entries.size() > firstNew ? TRUE : FALSE);
	if (mergedSourceContract) {
	    this->d->displayMeshLodContractRevisionValid = TRUE;
	    this->d->displayMeshLodContractSourceRevision =
		this->sourceRevision.getValue();
	    this->d->displayMeshLodContractInputsRevision =
		this->inputsRevision.getValue();
	}
	if (committed)
	    committed(context);
	SbBox3f bounds;
	if (!this->hasExactSourceBounds() &&
	    compact_index_source_bounds(index, bounds) && !bounds.isEmpty())
	    (void)this->setSourceBoundsState(TRUE, bounds.getMin(),
		bounds.getMax(), FALSE);
	if (!this->d->compactSelectedPaths.empty() && this->selected.getValue())
	    this->selected = FALSE;
    };
    try {
	for (const BObolCompactOccurrence &occurrence : occurrences) {
	    if (denyPartialMerge && changed == 1)
		throw std::bad_alloc();
	    if (!occurrence.geometry)
		continue;
	    const bool isOverview = BU_STR_EQUAL(
		occurrence.summary.recordRole.getString(), "lod-overview");
	    const int newTier = compact_geometry_tier(
		occurrence.summary.geometryKind.getString());
	    const char *p = database_source_skip_leading_slash(
		occurrence.summary.path.getString());
	    std::unordered_map<std::string, size_t>::iterator found =
		(p && p[0]) ? index.entryIndexByPath.find(p) :
		index.entryIndexByPath.end();
	    if (found != index.entryIndexByPath.end()) {
		const BObolCompactInstanceEntry &existing =
		    index.entries[found->second];
		const int oldTier = compact_geometry_tier(
		    existing.shapeSummary.geometryKind.getString());
		const bool evolvingOverview =
		BU_STR_EQUAL(occurrence.summary.recordRole.getString(),
		    "lod-overview") &&
		BU_STR_EQUAL(existing.shapeSummary.recordRole.getString(),
		    "lod-overview");
		const bool richerDataContract =
		newTier == oldTier &&
		((occurrence.sourceMeshRequestValid &&
		  !existing.sourceMeshRequestValid) ||
		 (occurrence.lodBacked && !existing.lodBacked) ||
		 evolvingOverview || authoritativeGeometry);
		if (newTier < oldTier ||
		(newTier == oldTier && !richerDataContract))
		    continue;
		const bool replaceGeometry = evolvingOverview ||
		authoritativeGeometry;
		const bool merged = newTier == oldTier && !replaceGeometry ?
		compact_index_merge_source_contract(index, found->second,
		    occurrence) :
		compact_index_replace_entry_geometry(this, index,
		    found->second, occurrence, *this->d);
		if (merged) {
		    overviewArrived = overviewArrived || isOverview;
		    changed++;
		    changedEntries.push_back(found->second);
		    if (occurrence.sourceMeshRequestValid)
		        mergedSourceContract = true;
		}
		continue;
	    }
	    compact_occurrence_build input;
	    input.occurrence = occurrence;
	    input.semantic = compact_semantic_from_summary(occurrence.summary);
	    int unsupported = 0;
	    const size_t before = index.entries.size();
	    compact_add_occurrence(this, index, input, ordinal, unsupported,
		this->d.get());
	    if (unsupported || index.entries.size() == before)
		continue;
	    overviewArrived = overviewArrived || isOverview;
	    changed++;
	    changedEntries.push_back(index.entries.size() - 1);
	    if (occurrence.sourceMeshRequestValid)
		mergedSourceContract = true;
	}
    } catch (...) {
	/* Every completed occurrence is internally consistent. Publish that
	 * prefix even when preparation of a later occurrence fails. */
	publish();
	throw;
    }
    if (batchCompleted)
	*batchCompleted = TRUE;
    publish();
    if (!changed)
	return 0;
    if (getenv("BOBOL_LOD_TRACE_SOURCE_CONTRACT"))
	bu_log("BObol LoD source contract merged compact occurrences path=%s "
	       "changed=%d entries=%zu requests=%zu resident_progressive=%zu "
	       "targets=%zu\n",
	       this->path.getValue().getString(), changed,
	       index.entries.size(), index.sourceMeshRequestCount,
	       index.residentProgressiveGeometryCount,
	       index.displayLodTargetCount);

    return changed;
}

SbBool
SoBRLDatabaseSource::acceptsCompactStreamCertification(
    size_t expectedCount,
    const BObolCompactSourceProfile *profile)
    const
{
    if (!expectedCount || (profile && !profile->isValid(expectedCount)) ||
	(this->d->compactExpectedInstanceCountCertified &&
	 this->d->compactExpectedInstanceCount != expectedCount) ||
	(this->d->compactSourceProfile.valid &&
	 !this->d->compactSourceProfile.isValid(expectedCount)))
	return FALSE;
    if (profile && this->d->compactSourceProfile.valid) {
	const BObolCompactSourceProfile &current =
	    this->d->compactSourceProfile;
	if (current.occurrenceCount != profile->occurrenceCount ||
	    current.uniqueAssetCount != profile->uniqueAssetCount ||
	    current.encodedSourceBytes != profile->encodedSourceBytes ||
	    current.largestAssetBytes != profile->largestAssetBytes ||
	    current.reusedOccurrenceCount != profile->reusedOccurrenceCount)
	    return FALSE;
    }
    return TRUE;
}

int
SoBRLDatabaseSource::certifyCompactStream(size_t expectedCount,
    const BObolCompactSourceProfile *profile)
{
    if (!this->acceptsCompactStreamCertification(expectedCount, profile))
	return -1;

    this->d->compactExpectedInstanceCount = expectedCount;
    this->d->compactExpectedInstanceCountCertified = TRUE;
    if (profile)
	this->d->compactSourceProfile = *profile;
    /* The owner can learn the certified count after it has already drained
     * the terminal leaf batch.  Schedule retirement here as well as in the
     * merge path so producer/consumer ordering cannot leave the overview
     * permanent. */
    if (this->d->compactIndex && this->d->compactOverviewState ==
	    BObolCompactOccurrenceRegistryState::OverviewState::Visible &&
	compact_complete_stream_leaf_frontier(*this->d->compactIndex,
	    this->d->compactExpectedInstanceCount))
	this->d->compactOverviewState =
	    BObolCompactOccurrenceRegistryState::OverviewState::RetirementPending;
    return 1;
}

int
SoBRLDatabaseSource::reserveCompactStreamCapacity(size_t expectedCount)
{
    if (!expectedCount)
	return 0;
    if (this->d->compactExpectedInstanceCountCertified &&
	this->d->compactExpectedInstanceCount != expectedCount)
	return -1;
    if (!this->d->compactIndex) {
	std::unique_ptr<BObolCompactInstanceIndex> index(
	    new BObolCompactInstanceIndex);
	index->instances.reserve(expectedCount);
	this->d->compactIndex = index.release();
	this->d->compactIndexActive = TRUE;
	this->d->compactOccurrenceRegistry = TRUE;
    }
    BObolCompactInstanceIndex &index = *this->d->compactIndex;
    /* The expected count is for this exact streamed draw target, not its
     * containing tree root. Reserving before the first leaf merge avoids a
     * later vector doubling and its temporary duplicate instance storage. */
    if (index.instances.capacity() < expectedCount)
	index.instances.reserve(expectedCount);
    return 1;
}


int
SoBRLDatabaseSource::adoptCompactOccurrencesFrom(
    const SoBRLDatabaseSource *const *sources, size_t count,
    BObolSourceAdoptionEffects *effects)
{
    if ((!sources && count) || !count)
	return 0;

    /* Retain every original owner through observer delivery. An observer may
     * remove or replace either side after the complete commit. */
    SbModernUtils::SoNodeRef targetOwner(this);
    std::vector<SbModernUtils::SoNodeRef> sourceOwners;
    sourceOwners.reserve(count);
    size_t occurrenceCount = 0;
    for (size_t i = 0; i < count; ++i) {
	const SoBRLDatabaseSource *source = sources[i];
	if (!source || source == this || !source->d->compactIndex)
	    continue;
	const size_t sourceCount = source->d->compactIndex->entries.size();
	if (sourceCount > std::numeric_limits<size_t>::max() - occurrenceCount)
	    throw std::length_error("compact adoption population overflow");
	sourceOwners.emplace_back(const_cast<SoBRLDatabaseSource *>(source));
	occurrenceCount += sourceCount;
    }

    SbMatrix placementInverse = SbMatrix::identity();
    if (this->drawMatrixValid.getValue())
	placementInverse = this->drawMatrix.getValue().inverse();

    std::vector<BObolCompactOccurrence> occurrences;
    occurrences.reserve(occurrenceCount);
    for (const auto &owner : sourceOwners) {
	const auto *source = static_cast<const SoBRLDatabaseSource *>(owner.get());
	for (const BObolCompactInstanceEntry &entry :
	     source->d->compactIndex->entries) {
	    if (!entry.geometry)
		continue;

	    BObolCompactOccurrence occurrence;
	    occurrence.geometry = entry.geometry;
	    occurrence.summary = entry.shapeSummary;
	    occurrence.summary.path = entry.semantic.path;
	    occurrence.summary.sourceName = entry.semantic.sourceName;
	    occurrence.summary.visible = entry.authoredVisible;
	    occurrence.summary.selected = entry.selected;
	    occurrence.summary.highlighted = entry.highlighted;
	    occurrence.geometryTransform = entry.geometryTransform;
	    occurrence.localTransform = entry.placementTransform;
	    if (source->drawMatrixValid.getValue())
		occurrence.localTransform.multRight(
		    source->drawMatrix.getValue());
	    occurrence.localTransform.multRight(placementInverse);
	    occurrence.viewDependentCsgGeometry =
		entry.viewDependentCsgGeometry;
	    occurrence.lodBacked = entry.lodBacked;
	    occurrence.sourceMeshRequestValid = entry.sourceMeshRequestValid;
	    occurrence.sourceMeshRequest = entry.sourceMeshRequest;
	    occurrence.occurrenceIndex = entry.occurrenceIndex;
	    occurrence.booleanOperation = entry.booleanOperation;
	    occurrences.push_back(std::move(occurrence));
	}
    }

    int changed = 0;
    std::unique_ptr<BObolCompactInstanceIndex> candidate;
    BObolCompactOccurrenceRegistryState::OverviewState nextOverview =
	this->d->compactOverviewState;
    SbBool nextContractValid = this->d->displayMeshLodContractRevisionValid;
    uint32_t nextContractSourceRevision =
	this->d->displayMeshLodContractSourceRevision;
    uint32_t nextContractInputsRevision =
	this->d->displayMeshLodContractInputsRevision;
    std::vector<size_t> changedEntries;
    bool sparseChangesValid = false;
    bool coverageInvalidated = false;
    if (!occurrences.empty()) {
	SoBRLDatabaseSource *detached = this->createDetachedRealizationTemplate();
	SbModernUtils::SoNodeRef detachedOwner(detached);
	detached->unref();
	auto &next = *detached->d;
	const auto &current = *this->d;
	next.compactHandleSourceId = current.compactHandleSourceId;
	next.compactExpectedInstanceCount = current.compactExpectedInstanceCount;
	next.compactExpectedInstanceCountCertified =
	    current.compactExpectedInstanceCountCertified;
	next.compactOverviewState = current.compactOverviewState;
	next.compactVisibilityFrontierActive =
	    current.compactVisibilityFrontierActive;
	next.compactVisibilityFrontierDefault =
	    current.compactVisibilityFrontierDefault;
	next.compactVisibilityFrontier = current.compactVisibilityFrontier;
	next.compactVisibilityFrontierStates =
	    current.compactVisibilityFrontierStates;
	next.compactSelectedPaths = current.compactSelectedPaths;
	next.compactPresentationOverrides = current.compactPresentationOverrides;
	next.displayMeshLodContractRevisionValid =
	    current.displayMeshLodContractRevisionValid;
	next.displayMeshLodContractSourceRevision =
	    current.displayMeshLodContractSourceRevision;
	next.displayMeshLodContractInputsRevision =
	    current.displayMeshLodContractInputsRevision;
	if (current.compactIndex) {
	    next.compactIndex =
		new BObolCompactInstanceIndex(*current.compactIndex);
	    next.compactIndexActive = current.compactIndexActive;
	    next.compactOccurrenceRegistry = current.compactOccurrenceRegistry;
	}
	const uint64_t journalRevision = next.cadBatchRevision;
	const size_t previousCount = next.compactIndex ?
	    next.compactIndex->entries.size() : 0;
	changed = detached->mergeCompactOccurrences(occurrences, FALSE);
	if (changed > 0) {
	    candidate.reset(next.compactIndex);
	    next.compactIndex = nullptr;
	    nextOverview = next.compactOverviewState;
	    nextContractValid = next.displayMeshLodContractRevisionValid;
	    nextContractSourceRevision =
		next.displayMeshLodContractSourceRevision;
	    nextContractInputsRevision =
		next.displayMeshLodContractInputsRevision;
	    sparseChangesValid = detached->getCadBatchChangedEntries(
		journalRevision, changedEntries) != FALSE;
	    coverageInvalidated = candidate &&
		candidate->entries.size() > previousCount;
	}
    }

    if (changed <= 0) {
	if (!effects)
	    return 0;
	effects->prepare();
	effects->commit();
	std::exception_ptr failure;
	effects->notify(failure);
	if (failure) std::rethrow_exception(failure);
	return 0;
    }

    SbBox3f candidateBounds;
    const bool updateBounds = !this->hasExactSourceBounds() && candidate &&
	compact_index_source_bounds(*candidate, candidateBounds) &&
	!candidateBounds.isEmpty();
    const SbBool nextBoundsValid = updateBounds ? TRUE :
	this->sourceBoundsValid.getValue();
    const SbBool nextBoundsExact = updateBounds ? FALSE :
	this->sourceBoundsExact.getValue();
    const SbVec3f nextBoundsMinimum = updateBounds ? candidateBounds.getMin() :
	this->sourceBoundsMin.getValue();
    const SbVec3f nextBoundsMaximum = updateBounds ? candidateBounds.getMax() :
	this->sourceBoundsMax.getValue();
    const bool clearSelection = !this->d->compactSelectedPaths.empty() &&
	this->selected.getValue();
    PreparedFieldNotifications<5> notifications(*this, {{
	{&this->sourceBoundsValid,
	 this->sourceBoundsValid.getValue() != nextBoundsValid},
	{&this->sourceBoundsExact,
	 this->sourceBoundsExact.getValue() != nextBoundsExact},
	{&this->sourceBoundsMin,
	 !database_source_vec3f_equal(this->sourceBoundsMin.getValue(),
	     nextBoundsMinimum)},
	{&this->sourceBoundsMax,
	 !database_source_vec3f_equal(this->sourceBoundsMax.getValue(),
	     nextBoundsMaximum)},
	{&this->selected, clearSelection}
    }});
    if (effects)
	effects->prepare();

    BObolCompactInstanceIndex *previous = this->d->compactIndex;
    this->d->compactIndex = candidate.release();
    this->d->compactIndexActive = TRUE;
    this->d->compactOccurrenceRegistry = TRUE;
    this->d->compactOverviewState = nextOverview;
    this->d->displayMeshLodContractRevisionValid = nextContractValid;
    this->d->displayMeshLodContractSourceRevision =
	nextContractSourceRevision;
    this->d->displayMeshLodContractInputsRevision =
	nextContractInputsRevision;
    this->markCompiledAssemblyDirty();
    if (sparseChangesValid) {
	this->markCadBatchDirty(changedEntries);
	this->markDisplayMeshLodDirty(changedEntries,
	    coverageInvalidated ? TRUE : FALSE);
    } else {
	this->markCadBatchDirty();
	this->markDisplayMeshLodDirty();
    }
    if (updateBounds)
	(void)this->setSourceBoundsState(TRUE, nextBoundsMinimum,
	    nextBoundsMaximum, FALSE);
    if (clearSelection)
	this->selected = FALSE;
    delete previous;
    if (effects)
	effects->commit();

    notifications.restore();
    std::exception_ptr failure;
    notifications.notify(failure);
    if (effects)
	effects->notify(failure);
    if (failure) std::rethrow_exception(failure);
    return changed;
}

static bool
database_source_retarget_path_component(std::string &path,
	const std::string &oldComponent, const std::string &newComponent)
{
    if (path.empty() || oldComponent.empty() || newComponent.empty())
	return false;

    bool changed = false;
    size_t start = 0;
    while (start < path.size()) {
	const size_t end = path.find('/', start);
	const size_t componentEnd =
	    end == std::string::npos ? path.size() : end;
	const size_t suffix = path.find('@', start);
	const size_t baseEnd =
	    suffix != std::string::npos && suffix < componentEnd ?
	    suffix : componentEnd;
	if (path.compare(start, baseEnd - start, oldComponent) == 0 &&
	    baseEnd - start == oldComponent.size()) {
	    path.replace(start, oldComponent.size(), newComponent);
	    const size_t nextEnd =
		path.find('/', start + newComponent.size());
	    start = nextEnd == std::string::npos ? path.size() : nextEnd;
	    changed = true;
	} else {
	    start = componentEnd;
	}
	if (start < path.size() && path[start] == '/')
	    start++;
    }
    return changed;
}

bool
bobol_database_retarget_path_components(
    std::string &path, const char *oldObjectPath, const char *newObjectPath)
{
    if (!oldObjectPath || !oldObjectPath[0] ||
	!newObjectPath || !newObjectPath[0])
	return false;

    const std::string oldNormalized =
	database_source_skip_leading_slash(oldObjectPath);
    const std::string newNormalized =
	database_source_skip_leading_slash(newObjectPath);
    const size_t oldSlash = oldNormalized.find_last_of('/');
    const size_t newSlash = newNormalized.find_last_of('/');
    const std::string oldComponent = oldNormalized.substr(
	oldSlash == std::string::npos ? 0 : oldSlash + 1);
    const std::string newComponent = newNormalized.substr(
	newSlash == std::string::npos ? 0 : newSlash + 1);
    if (oldComponent.empty() || newComponent.empty() ||
	oldComponent == newComponent)
	return false;
    return database_source_retarget_path_component(path, oldComponent,
	newComponent);
}

static int
compact_retarget_occurrence_paths(BObolCompactInstanceIndex &index,
    const char *oldPath, const char *newPath)
{
    const std::string oldNormalized =
	database_source_skip_leading_slash(oldPath);
    const std::string newNormalized =
	database_source_skip_leading_slash(newPath);
    const size_t oldSlash = oldNormalized.find_last_of('/');
    const size_t newSlash = newNormalized.find_last_of('/');
    const std::string oldComponent = oldNormalized.substr(
	oldSlash == std::string::npos ? 0 : oldSlash + 1);
    const std::string newComponent = newNormalized.substr(
	newSlash == std::string::npos ? 0 : newSlash + 1);
    if (oldComponent.empty() || newComponent.empty() ||
	oldComponent == newComponent)
	return 0;

    int changed = 0;
    for (size_t i = 0; i < index.entries.size(); i++) {
	BObolCompactInstanceEntry &entry = index.entries[i];
	std::string semanticPath = entry.semantic.path.getString();
	if (!database_source_retarget_path_component(semanticPath,
		oldComponent, newComponent))
	    continue;

	entry.semantic.path = semanticPath.c_str();
	std::string sourceName = entry.semantic.sourceName.getString();
	if (database_source_retarget_path_component(sourceName,
		oldComponent, newComponent))
	    entry.semantic.sourceName = sourceName.c_str();
	compact_note_semantic_change(entry);
	if (entry.sourceMeshRequestValid) {
	    std::string requestPath =
		entry.sourceMeshRequest.path.getString();
	    (void)database_source_retarget_path_component(requestPath,
		oldComponent, newComponent);
	    entry.sourceMeshRequest.path = requestPath.c_str();
	    std::string assetPath =
		entry.sourceMeshRequest.meshAssetPath.getString();
	    (void)database_source_retarget_path_component(assetPath,
		oldComponent, newComponent);
	    entry.sourceMeshRequest.meshAssetPath = assetPath.c_str();
	    std::string requestName =
		entry.sourceMeshRequest.sourceName.getString();
	    if (database_source_retarget_path_component(requestName,
		    oldComponent, newComponent))
		entry.sourceMeshRequest.sourceName = requestName.c_str();
	    std::string assetName =
		entry.sourceMeshRequest.meshAssetName.getString();
	    if (database_source_retarget_path_component(assetName,
		    oldComponent, newComponent))
		entry.sourceMeshRequest.meshAssetName = assetName.c_str();
	}
	if (i < index.instances.size())
	    index.instances[i].record.childName = semanticPath;
	compact_sync_shape_summary_state(entry);
	changed++;
    }
    if (changed)
	compact_rebuild_entry_index(index);
    return changed;
}


int
SoBRLDatabaseSource::retargetCompactOccurrencePaths(
    const char *oldPath, const char *newPath)
{
    if (!this->d->compactIndex || !oldPath || !oldPath[0] ||
	!newPath || !newPath[0])
	return 0;
    BObolSourceRenamePublication publication(*this, oldPath, newPath);
    if (!publication.changed())
	return 0;
    publication.commit();
    publication.restore();
    std::exception_ptr failure;
    publication.notify(failure);
    if (failure)
	std::rethrow_exception(failure);
    return publication.compactPathChangeCount();
}

static uint64_t
compact_structure_signature(const BObolCompactInstanceIndex *index)
{
    if (!index)
	return 0;
    uint64_t signature = 1469598103934665603ULL;
    auto mix = [&signature](uint64_t value) {
	signature ^= value;
	signature *= 1099511628211ULL;
    };
    mix(index->parts.size());
    mix(index->entries.size());
    for (const BObolCompactInstanceEntry &entry : index->entries) {
	mix(entry.instance.w0);
	mix(entry.instance.w1);
	mix(entry.part.w0);
	mix(entry.part.w1);
	mix(entry.geometryRevision);
	mix(entry.placementRevision);
    }
    return signature ? signature : 1;
}

static void
compact_signature_mix(uint64_t &signature, uint64_t value)
{
    signature ^= value;
    signature *= 1099511628211ULL;
}

static void
compact_signature_mix_string(uint64_t &signature, const SbString &value)
{
    const char *string = value.getString();
    const size_t length = string ? strlen(string) : 0;
    compact_signature_mix(signature, length);
    for (size_t i = 0; i < length; i++)
	compact_signature_mix(signature, static_cast<unsigned char>(string[i]));
}

static void
compact_signature_mix_float(uint64_t &signature, float value)
{
    uint32_t bits = 0;
    memcpy(&bits, &value, sizeof(bits));
    compact_signature_mix(signature, bits);
}

/* Pick metadata is independent of the retained geometry and visual state.
 * Updating the semantic map rewrites a sorted map and copies several strings
 * for every occurrence, so keep it out of style/visibility-only updates. */
static uint64_t
compact_semantic_signature(const BObolCompactInstanceIndex *index)
{
    if (!index)
	return 0;

    uint64_t signature = 1469598103934665603ULL;
    compact_signature_mix(signature, index->entries.size());
    for (const BObolCompactInstanceEntry &entry : index->entries) {
	const SoBRLCadAssembly::InstanceSemantic &semantic = entry.semantic;
	compact_signature_mix(signature, entry.instance.w0);
	compact_signature_mix(signature, entry.instance.w1);
	/* setInstanceSemantic publishes this entry value, not the semantic's
	 * stale construction-time value. */
	compact_signature_mix_string(signature, entry.instanceKey);
	compact_signature_mix_string(signature, semantic.path);
	compact_signature_mix_string(signature, semantic.sourceName);
	compact_signature_mix_string(signature, semantic.sourceType);
	compact_signature_mix_string(signature, semantic.materialShader);
	compact_signature_mix_string(signature, semantic.editIntentId);
	compact_signature_mix_string(signature, semantic.editIntentRole);
	compact_signature_mix(signature, semantic.sourceId);
	compact_signature_mix(signature,
	    static_cast<uint32_t>(semantic.regionId));
	compact_signature_mix(signature,
	    static_cast<uint32_t>(semantic.airCode));
	compact_signature_mix(signature,
	    static_cast<uint32_t>(semantic.materialId));
	compact_signature_mix(signature, static_cast<uint32_t>(semantic.los));
	compact_signature_mix(signature, semantic.materialColorValid ? 1 : 0);
	if (semantic.materialColorValid) {
	    compact_signature_mix_float(signature, semantic.materialColor[0]);
	    compact_signature_mix_float(signature, semantic.materialColor[1]);
	    compact_signature_mix_float(signature, semantic.materialColor[2]);
	}
	compact_signature_mix(signature,
	    static_cast<uint32_t>(semantic.primitiveKind));
    }
    return signature ? signature : 1;
}

/* The retained assembly owns a copy of every instance style.  Track the
 * compact revisions which can change that copy so an unrelated source-node
 * notification (for example, a camera redraw) does not republish every style
 * and touch the retained assembly again. */
static uint64_t
compact_style_signature(const BObolCompactInstanceIndex *index)
{
    if (!index)
	return 0;

    uint64_t signature = 1469598103934665603ULL;
    compact_signature_mix(signature, index->entries.size());
    for (const BObolCompactInstanceEntry &entry : index->entries) {
	compact_signature_mix(signature, entry.instance.w0);
	compact_signature_mix(signature, entry.instance.w1);
	compact_signature_mix(signature, entry.appearanceRevision);
	compact_signature_mix(signature, entry.visibilityRevision);
	compact_signature_mix(signature, entry.selectionRevision);
    }
    return signature ? signature : 1;
}

static void
compact_assembly_draw_mode(SoBRLCadAssembly *assembly,
    const SoBRLDatabaseSource *source,
    const BObolCompactInstanceIndex *index)
{
    if (!assembly || !source || !index)
	return;
    assembly->setPresentationDrawMode(cad_presentation_draw_mode(
	source_record_draw_mode(source), index->shadedCount,
	index->wireCount));
}

uint64_t
SoBRLDatabaseSource::cadBatchStructureSignature(void) const
{
    uint64_t signature = this->d->compactIndexActive && this->d->compactIndex ?
	compact_structure_signature(this->d->compactIndex) : this->d->cadBatchRevision;
    signature ^= static_cast<uint64_t>(source_record_draw_mode(this));
    signature *= 1099511628211ULL;
    signature ^= static_cast<uint64_t>(this->visible.getValue());
    return signature ? signature : 1;
}

uint64_t
SoBRLDatabaseSource::cadBatchStyleSignature(void) const
{
    return this->d->compactIndexActive && this->d->compactIndex ?
	compact_style_signature(this->d->compactIndex) :
	this->d->cadBatchRevision;
}

uint64_t
SoBRLDatabaseSource::cadBatchSemanticSignature(void) const
{
    return this->d->compactIndexActive && this->d->compactIndex ?
	compact_semantic_signature(this->d->compactIndex) : 0;
}

static std::vector<BObolCadPresentationBridgeState::CompactStamp>
compact_presentation_stamps(const BObolCompactInstanceIndex *index)
{
    std::vector<BObolCadPresentationBridgeState::CompactStamp> stamps;
    if (!index)
	return stamps;
    stamps.reserve(index->entries.size());
    for (const BObolCompactInstanceEntry &entry : index->entries) {
	BObolCadPresentationBridgeState::CompactStamp stamp;
	stamp.instance = entry.instance;
	stamp.part = entry.part;
	stamp.geometryRevision = entry.geometryRevision;
	stamp.placementRevision = entry.placementRevision;
	stamp.appearanceRevision = entry.appearanceRevision;
	stamp.visibilityRevision = entry.visibilityRevision;
	stamp.selectionRevision = entry.selectionRevision;
	stamp.semanticRevision = entry.semanticRevision;
	stamps.push_back(stamp);
    }
    return stamps;
}

static bool
compact_structure_stamps_equal(
    const std::vector<BObolCadPresentationBridgeState::CompactStamp> &a,
    const std::vector<BObolCadPresentationBridgeState::CompactStamp> &b)
{
    if (a.size() != b.size())
	return false;
    for (size_t i = 0; i < a.size(); ++i) {
	if (!a[i].sameStructure(b[i]))
	    return false;
    }
    return true;
}

static bool
compact_style_stamps_equal(
    const std::vector<BObolCadPresentationBridgeState::CompactStamp> &a,
    const std::vector<BObolCadPresentationBridgeState::CompactStamp> &b)
{
    if (a.size() != b.size())
	return false;
    for (size_t i = 0; i < a.size(); ++i) {
	if (!a[i].sameStyle(b[i]))
	    return false;
    }
    return true;
}

static bool
compact_semantic_stamps_equal(
    const std::vector<BObolCadPresentationBridgeState::CompactStamp> &a,
    const std::vector<BObolCadPresentationBridgeState::CompactStamp> &b)
{
    if (a.size() != b.size())
	return false;
    for (size_t i = 0; i < a.size(); ++i) {
	if (!a[i].sameSemantic(b[i]))
	    return false;
    }
    return true;
}

int
SoBRLDatabaseSource::syncCompiledAssembly(void)
{
    const SbUniqueId sourceNodeId = this->getNodeId();
    if (!this->d->compiledAssemblyDirty &&
	this->d->compiledAssemblyNodeId == sourceNodeId)
	return this->d->compiledAssemblyActive ? 1 : 0;
    const SbBool precedingAssemblyActive = this->d->compiledAssemblyActive;
    const auto rejectPublication = [this, precedingAssemblyActive]() {
	this->d->compiledAssemblyActive = precedingAssemblyActive;
	return precedingAssemblyActive ? 1 : 0;
    };
    this->d->compiledAssemblyActive = FALSE;

    const SbBool hasCompactPayload = this->d->compactIndexActive &&
	this->d->compactIndex && !this->d->compactIndex->instances.empty();
    if (!this->visible.getValue() ||
	this->auxiliarySource.getValue() ||
	(!hasCompactPayload &&
	 (this->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	  this->needsRealization())) ||
	source_has_auxiliary_children(this)) {
	this->d->clearCompiledCompactEvidence();
	this->d->compiledAssemblyNodeId = sourceNodeId;
	this->d->compiledAssemblyDirty = FALSE;
	return 0;
    }

    if (!this->d->compiledAssembly) {
	this->d->compiledAssembly = new SoBRLCadAssembly;
	this->d->compiledAssembly->ref();
    }

    const std::vector<BObolCadPresentationBridgeState::CompactStamp>
	compactStamps = hasCompactPayload ?
	    compact_presentation_stamps(this->d->compactIndex) :
	    std::vector<BObolCadPresentationBridgeState::CompactStamp>();
    const bool structureCurrent = hasCompactPayload &&
	this->d->compiledCompactPartCount ==
	    this->d->compactIndex->parts.size() &&
	compact_structure_stamps_equal(
	    this->d->compiledCompactStamps, compactStamps) &&
	this->d->compiledAssembly->instanceCount() ==
	    this->d->compactIndex->instances.size() &&
	this->d->compiledAssembly->partCount() ==
	    this->d->compactIndex->parts.size();
    if (structureCurrent) {
	if (!compact_style_stamps_equal(
		this->d->compiledCompactStamps, compactStamps)) {
	    std::vector<Obol::InstanceStyleUpdate> styles;
	    styles.reserve(this->d->compactIndex->instances.size());
	for (const Obol::InstanceUpdate &instance :
		 this->d->compactIndex->instances) {
		Obol::InstanceStyleUpdate style;
		style.instance = instance.instance;
		style.style = instance.record.style;
	    styles.push_back(style);
	}
	    if (!bobol_cad_validate_styles(styles,
		    "compiled CAD style preflight") ||
		!bobol_cad_publish_styles(this->d->compiledAssembly, styles,
		    "compiled CAD style publication"))
		return rejectPublication();
	}
	if (this->d->compiledCompactHidden !=
		this->d->compactIndex->hiddenInstances) {
	    this->d->compiledAssembly->setHiddenInstances(
		this->d->compactIndex->hiddenInstances);
	}
	if (this->d->compiledCompactSelected !=
		this->d->compactIndex->selectedInstances) {
	    this->d->compiledAssembly->setSelectedInstances(
		this->d->compactIndex->selectedInstances);
	}
	if (this->d->compiledCompactUnpickable !=
		this->d->compactIndex->unpickableInstances) {
	    this->d->compiledAssembly->setUnpickableInstances(
		this->d->compactIndex->unpickableInstances);
	}
	if (!compact_semantic_stamps_equal(
		this->d->compiledCompactStamps, compactStamps)) {
	    for (const BObolCompactInstanceEntry &entry :
		 this->d->compactIndex->entries) {
		SoBRLCadAssembly::InstanceSemantic semantic = entry.semantic;
		semantic.sourceInstanceKey = compact_instance_identity(entry);
		this->d->compiledAssembly->setInstanceSemantic(entry.instance, semantic);
	    }
	}
	this->d->compiledCompactStamps = compactStamps;
	this->d->compiledCompactHidden =
	    this->d->compactIndex->hiddenInstances;
	this->d->compiledCompactSelected =
	    this->d->compactIndex->selectedInstances;
	this->d->compiledCompactUnpickable =
	    this->d->compactIndex->unpickableInstances;
	compact_assembly_draw_mode(this->d->compiledAssembly, this,
	    this->d->compactIndex);
	this->d->compiledAssemblyActive = TRUE;
	this->ensureCompiledAssemblyChild();
	this->d->compiledAssemblyNodeId = this->getNodeId();
	this->d->compiledAssemblyDirty = FALSE;
	return 1;
    }

    if (this->d->compactIndexActive && this->d->compactIndex &&
	!this->d->compactIndex->instances.empty()) {
	std::vector<Obol::PartUpdate> parts;
	parts.reserve(this->d->compactIndex->parts.size());
	for (const BObolCompactPartReference &partRef :
	     this->d->compactIndex->parts) {
	    if (!partRef.geometry)
		continue;
	    Obol::PartUpdate part;
	    part.part = partRef.part;
	    if (!bobol_cad_admit_geometry(partRef.geometry, part.geometry,
		    "compiled compact part staging"))
		return rejectPublication();
	    parts.push_back(part);
	}
	if (!bobol_cad_replace_scene(this->d->compiledAssembly, parts,
		this->d->compactIndex->instances,
		"compiled compact replacement"))
	    return rejectPublication();
	this->d->compiledAssembly->clearSemanticMap();
	this->d->compiledAssembly->setHiddenInstances(
	    this->d->compactIndex->hiddenInstances);
	this->d->compiledAssembly->setSelectedInstances(
	    this->d->compactIndex->selectedInstances);
	this->d->compiledAssembly->setUnpickableInstances(
	    this->d->compactIndex->unpickableInstances);
	for (size_t i = 0; i < this->d->compactIndex->entries.size(); i++) {
	    const BObolCompactInstanceEntry &entry =
		this->d->compactIndex->entries[i];
	    SoBRLCadAssembly::InstanceSemantic semantic = entry.semantic;
	    semantic.sourceInstanceKey = compact_instance_identity(entry);
	    this->d->compiledAssembly->setInstanceSemantic(entry.instance, semantic);
	}
	compact_assembly_draw_mode(this->d->compiledAssembly, this,
	    this->d->compactIndex);
	this->d->compiledAssemblyActive = TRUE;
	this->d->compiledCompactStamps = compactStamps;
	this->d->compiledCompactPartCount =
	    this->d->compactIndex->parts.size();
	this->d->compiledCompactHidden =
	    this->d->compactIndex->hiddenInstances;
	this->d->compiledCompactSelected =
	    this->d->compactIndex->selectedInstances;
	this->d->compiledCompactUnpickable =
	    this->d->compactIndex->unpickableInstances;
	this->ensureCompiledAssemblyChild();
	this->d->compiledAssemblyNodeId = this->getNodeId();
	this->d->compiledAssemblyDirty = FALSE;
	return 1;
    }

    /* The retained compact signature has no meaning for a child-graph build. */
    this->d->clearCompiledCompactEvidence();
    cad_build_data data;
    data.source = this;
    data.ordinal = 0;
    data.unsupported = 0;
    data.wireCount = 0;
    data.shadedCount = 0;

    const SbMatrix identity = SbMatrix::identity();
    for (int i = 0; i < this->getNumChildren() && !data.unsupported; i++)
	cad_collect_realized_node(data, this->getChild(i), identity);

    std::vector<Obol::PartUpdate> parts;
    if (!data.unsupported && !data.instances.empty()) {
	parts.reserve(data.parts.size());
	for (cad_pending_part &partUpdate : data.parts) {
	    Obol::PartUpdate part;
	    part.part = partUpdate.part;
	    const auto geometry =
		bobol_cad_build_geometry(std::move(partUpdate.geometry),
		    "compiled realized part staging");
	    if (!bobol_cad_admit_geometry(geometry, part.geometry,
		    "compiled realized part staging"))
		return rejectPublication();
	    parts.push_back(std::move(part));
	}
    }

    const std::vector<Obol::InstanceUpdate> emptyInstances;
    const std::vector<Obol::InstanceUpdate> &replacementInstances =
	(!data.unsupported && !data.instances.empty()) ?
	data.instances : emptyInstances;
    if (!bobol_cad_replace_scene(this->d->compiledAssembly, parts,
	    replacementInstances, "compiled realized replacement"))
	return rejectPublication();
    this->d->compiledAssembly->clearSemanticMap();
    if (!data.unsupported && !data.instances.empty()) {
	for (const auto &semantic : data.semantics)
	    this->d->compiledAssembly->setInstanceSemantic(
		semantic.first, semantic.second);
	this->d->compiledAssembly->setHiddenInstances(data.hiddenInstances);
	this->d->compiledAssembly->setSelectedInstances(data.selectedInstances);
	this->d->compiledAssembly->setUnpickableInstances(
	    data.unpickableInstances);
	this->d->compiledAssembly->setPresentationDrawMode(
	    cad_presentation_draw_mode(source_record_draw_mode(this),
		data.shadedCount, data.wireCount));
	this->d->compiledAssemblyActive = TRUE;
    }
    this->d->compiledAssemblyNodeId = this->getNodeId();
    this->d->compiledAssemblyDirty = FALSE;
    return this->d->compiledAssemblyActive ? 1 : 0;
}

int
SoBRLDatabaseSource::appendCadRenderBatch(BObolCadBatchBuildState *state,
	SbBool includeGeometry, SbBool includeSemantics)
{
    const SbBool hasCompactPayload = this->d->compactIndexActive &&
	this->d->compactIndex && !this->d->compactIndex->instances.empty();
    if (!state || !state->valid || !this->visible.getValue() ||
	this->auxiliarySource.getValue() ||
	source_record_draw_mode(this) == BOBOL_LOD_DRAW_HIDDEN_LINE ||
	(!hasCompactPayload &&
	 (this->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	  this->needsRealization())) || source_has_auxiliary_children(this))
	return 0;

    if (this->d->compactIndexActive && this->d->compactIndex &&
	!this->d->compactIndex->instances.empty()) {
	if (includeGeometry) {
	    for (const BObolCompactPartReference &partRef :
		 this->d->compactIndex->parts) {
		if (!partRef.geometry ||
		    !state->partIds.insert(partRef.part).second)
		    continue;
		Obol::PartUpdate part;
		part.part = partRef.part;
		if (!bobol_cad_admit_geometry(partRef.geometry, part.geometry,
			"CAD render batch staging")) {
		    state->valid = false;
		    return 0;
		}
		state->parts.push_back(part);
	    }
	}
	state->instances.insert(state->instances.end(),
	    this->d->compactIndex->instances.begin(),
	    this->d->compactIndex->instances.end());
	state->hiddenInstances.insert(state->hiddenInstances.end(),
	    this->d->compactIndex->hiddenInstances.begin(),
	    this->d->compactIndex->hiddenInstances.end());
	state->selectedInstances.insert(state->selectedInstances.end(),
	    this->d->compactIndex->selectedInstances.begin(),
	    this->d->compactIndex->selectedInstances.end());
	state->unpickableInstances.insert(state->unpickableInstances.end(),
	    this->d->compactIndex->unpickableInstances.begin(),
	    this->d->compactIndex->unpickableInstances.end());
	if (includeSemantics) {
	    for (const BObolCompactInstanceEntry &entry :
		 this->d->compactIndex->entries) {
		SoBRLCadAssembly::InstanceSemantic semantic = entry.semantic;
		semantic.sourceInstanceKey = compact_instance_identity(entry);
		state->semantics.emplace_back(entry.instance, semantic);
	    }
	}
	state->wireCount += this->d->compactIndex->wireCount;
	state->shadedCount += this->d->compactIndex->shadedCount;
	return 1;
    }

    cad_build_data data;
    data.source = this;
    data.ordinal = 0;
    data.unsupported = 0;
    data.wireCount = 0;
    data.shadedCount = 0;

    const SbMatrix identity = SbMatrix::identity();
    for (int i = 0; i < this->getNumChildren() && !data.unsupported; i++)
	cad_collect_realized_node(data, this->getChild(i), identity);
    if (data.unsupported || data.instances.empty())
	return 0;

    for (cad_pending_part &part : data.parts) {
	if (state->partIds.insert(part.part).second) {
	    Obol::PartUpdate sharedPart;
	    sharedPart.part = part.part;
	    const auto sharedGeometry =
		bobol_cad_build_geometry(std::move(part.geometry),
		    "collected CAD render batch staging");
	    if (!bobol_cad_admit_geometry(sharedGeometry, sharedPart.geometry,
		    "collected CAD render batch staging")) {
		state->valid = false;
		return 0;
	    }
	    state->parts.push_back(std::move(sharedPart));
	}
    }
    state->instances.insert(state->instances.end(), data.instances.begin(),
	data.instances.end());
    if (includeSemantics)
	state->semantics.insert(state->semantics.end(), data.semantics.begin(),
	    data.semantics.end());
    state->hiddenInstances.insert(state->hiddenInstances.end(),
	data.hiddenInstances.begin(), data.hiddenInstances.end());
    state->selectedInstances.insert(state->selectedInstances.end(),
	data.selectedInstances.begin(), data.selectedInstances.end());
    state->unpickableInstances.insert(state->unpickableInstances.end(),
	data.unpickableInstances.begin(), data.unpickableInstances.end());
    state->wireCount += data.wireCount;
    state->shadedCount += data.shadedCount;
    return 1;
}

SoBRLDatabaseSource::SoBRLDatabaseSource(void) :
    SoBRLDatabaseSource(FieldObservation::Observe)
{}

SoBRLDatabaseSource::SoBRLDatabaseSource(FieldObservation observation) :
    d(new Impl)
{
    SO_NODE_CONSTRUCTOR(SoBRLDatabaseSource);

    SO_NODE_DEFINE_ENUM_VALUE(DrawMode, WIREFRAME);
    SO_NODE_DEFINE_ENUM_VALUE(DrawMode, SHADED);
    SO_NODE_DEFINE_ENUM_VALUE(RealizationStatus, UNREALIZED);
    SO_NODE_DEFINE_ENUM_VALUE(RealizationStatus, REALIZED);
    SO_NODE_DEFINE_ENUM_VALUE(RealizationStatus, FAILED);
    SO_NODE_DEFINE_ENUM_VALUE(MaterialPolicy, MATERIAL_INHERIT);
    SO_NODE_DEFINE_ENUM_VALUE(MaterialPolicy, MATERIAL_DATABASE);

    SO_NODE_ADD_FIELD(instanceKey, (""));
    SO_NODE_ADD_FIELD(path, (""));
    SO_NODE_ADD_FIELD(parentInstanceKey, (""));
    SO_NODE_ADD_FIELD(occurrenceIndex, (0));
    SO_NODE_ADD_FIELD(booleanOperation, (BOOLEAN_UNION));
    SO_NODE_ADD_FIELD(displayName, (""));
    SO_NODE_ADD_FIELD(representationKey, (""));
    SO_NODE_ADD_FIELD(representationMode, (-1));
    SO_NODE_ADD_FIELD(auxiliarySource, (FALSE));
    SO_NODE_ADD_FIELD(drawMode, (WIREFRAME));
    SO_NODE_SET_SF_ENUM_TYPE(drawMode, DrawMode);
    SO_NODE_ADD_FIELD(visible, (TRUE));
    SO_NODE_ADD_FIELD(selected, (FALSE));
    SO_NODE_ADD_FIELD(highlighted, (FALSE));
    SO_NODE_ADD_FIELD(lineStyle, (0));
    SO_NODE_ADD_FIELD(lineWidth, (0));
    SO_NODE_ADD_FIELD(transparency, (0.0f));
    SO_NODE_ADD_FIELD(materialColorValid, (FALSE));
    SO_NODE_ADD_FIELD(materialColor, (SbColor(1.0f, 1.0f, 1.0f)));
    SO_NODE_ADD_FIELD(materialRevision, (0));
    SO_NODE_ADD_FIELD(materialPolicy, (MATERIAL_INHERIT));
    SO_NODE_SET_SF_ENUM_TYPE(materialPolicy, MaterialPolicy);
    SO_NODE_ADD_FIELD(databaseMetadataValid, (FALSE));
    SO_NODE_ADD_FIELD(databaseRegionId, (0));
    SO_NODE_ADD_FIELD(databaseAirCode, (0));
    SO_NODE_ADD_FIELD(databaseMaterialId, (0));
    SO_NODE_ADD_FIELD(databaseLos, (0));
    SO_NODE_ADD_FIELD(databaseMaterialColorValid, (FALSE));
    SO_NODE_ADD_FIELD(databaseMaterialColor,
		      (SbColor(1.0f, 1.0f, 1.0f)));
    SO_NODE_ADD_FIELD(databaseMaterialShader, (""));
    SO_NODE_ADD_FIELD(colorOverride, (FALSE));
    SO_NODE_ADD_FIELD(color, (SbColor(1.0f, 1.0f, 1.0f)));
    SO_NODE_ADD_FIELD(selectedColor, (SbColor(1.0f, 1.0f, 1.0f)));
    SO_NODE_ADD_FIELD(highlightedColor, (SbColor(1.0f, 1.0f, 0.0f)));
    SO_NODE_ADD_FIELD(ghostedColor, (SbColor(0.55f, 0.55f, 0.55f)));
    SO_NODE_ADD_FIELD(drawMatrixValid, (FALSE));
    SO_NODE_ADD_FIELD(drawMatrix, (SbMatrix::identity()));
    SO_NODE_ADD_FIELD(drawCenterValid, (FALSE));
    SO_NODE_ADD_FIELD(drawCenter, (SbVec3f(0.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_FIELD(drawSizeValid, (FALSE));
    SO_NODE_ADD_FIELD(drawSize, (0.0f));
    SO_NODE_ADD_FIELD(sourceBoundsValid, (FALSE));
    SO_NODE_ADD_FIELD(sourceBoundsExact, (FALSE));
    SO_NODE_ADD_FIELD(sourceBoundsMin, (SbVec3f(0.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_FIELD(sourceBoundsMax, (SbVec3f(0.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_FIELD(tessellationAbsTol, (0.0f));
    SO_NODE_ADD_FIELD(tessellationRelTol, (0.01f));
    SO_NODE_ADD_FIELD(tessellationNormTol, (0.0f));
    SO_NODE_ADD_FIELD(lodBotThreshold, (0));
    SO_NODE_ADD_FIELD(sourceRevision, (0));
    SO_NODE_ADD_FIELD(inputsRevision, (0));
    SO_NODE_ADD_FIELD(viewRevision, (0));
    SO_NODE_ADD_FIELD(realizedRevision, (0));
    SO_NODE_ADD_FIELD(realizedSourceRevision, (0));
    SO_NODE_ADD_FIELD(realizedInputsRevision, (0));
    SO_NODE_ADD_FIELD(realizedViewRevision, (0));
    SO_NODE_ADD_FIELD(realizationStatus, (UNREALIZED));
    SO_NODE_SET_SF_ENUM_TYPE(realizationStatus, RealizationStatus);
    SO_NODE_ADD_FIELD(realizationDiagnostic, (""));
    SO_NODE_ADD_FIELD(realizationIdentity, (""));
    SO_NODE_ADD_FIELD(realizationRoleFlags, (REALIZATION_ROLE_NONE));
    SO_NODE_ADD_FIELD(realizationViewDependent, (FALSE));
    SO_NODE_ADD_FIELD(realizationCsgLodEnabled, (FALSE));
    SO_NODE_ADD_FIELD(realizationMeshLodEnabled, (FALSE));
    SO_NODE_ADD_FIELD(realizationViewScale, (0.0f));
    SO_NODE_ADD_FIELD(realizationLodScale, (1.0f));
    SO_NODE_ADD_FIELD(realizationViewWidth, (0));
    SO_NODE_ADD_FIELD(realizationViewHeight, (0));
    SO_NODE_ADD_FIELD(realizationBotThreshold, (0));
    SO_NODE_ADD_FIELD(realizationCurveScale, (0.0f));
    SO_NODE_ADD_FIELD(realizationPointScale, (0.0f));
    SO_NODE_ADD_FIELD(stale, (TRUE));
    SO_NODE_ADD_FIELD(staleReason, (STALE_SOURCE));

    if (observation == FieldObservation::Observe) {
	this->d->observeFields = true;
	try { this->attachObservedFieldAuditors(); }
	catch (...) {
	    this->detachObservedFieldAuditors();
	    throw;
	}
    }
}

SoBRLDatabaseSource::~SoBRLDatabaseSource(void)
{
    this->d->observeFields = false;
    this->detachObservedFieldAuditors();
    this->clearCompiledAssembly();
    this->discardCompactInstanceHistory();
    this->clearMeshLod();
}

void
SoBRLDatabaseSource::initClass(void)
{
    SoCADAssembly::initClass();
    SoBRLCadAssembly::initClass();
    SoBRLCadRenderBatch::initClass();
    SO_NODE_INIT_CLASS(SoBRLDatabaseSource, SoSeparator, "Separator");
}

struct ObservedFieldChange {
    bool watched = false;
    uint32_t reason = SoBRLDatabaseSource::STALE_NONE;
};

static bool
source_external_realization_supports_draw_mode(
    const SoBRLDatabaseSource &source, int mode)
{
    return (source.realizationRoleFlags.getValue() &
	    SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL) &&
	source.realizationStatus.getValue() == SoBRLDatabaseSource::REALIZED &&
	!source.stale.getValue() &&
	((mode == SoBRLDatabaseSource::SHADED &&
	  source.hasRealizedMeshGeometry()) ||
	 (mode == SoBRLDatabaseSource::WIREFRAME &&
	  source.hasRealizedWireGeometry()));
}

static ObservedFieldChange
observed_field_change(const SoBRLDatabaseSource &source, const SoField *field)
{
    ObservedFieldChange change;
    change.watched = field == &source.path || field == &source.instanceKey ||
	field == &source.representationKey || field == &source.representationMode ||
	field == &source.drawMode || field == &source.tessellationAbsTol ||
	field == &source.tessellationRelTol || field == &source.tessellationNormTol ||
	field == &source.lodBotThreshold || field == &source.sourceRevision ||
	field == &source.inputsRevision || field == &source.viewRevision;
    if (!change.watched)
	return change;

    const int terminalEvaluated =
	source.representationMode.getValue() == SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE ||
	source.representationMode.getValue() == SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS;
    /* Evaluated sources own a complete terminal CSG result.  The camera and
     * draw-channel fields which drive progressive source production do not
     * change that result; configuration notifications for those fields must
	 * not send it back through the view-dependent pipeline. */
    if (terminalEvaluated &&
	(field == &source.viewRevision || field == &source.lodBotThreshold ||
	 field == &source.drawMode || field == &source.representationMode))
	return change;

    /* A delayed revision callback for an already realized epoch adds no
     * invalidation; configuration changes still synchronize their owners. */
    change.reason = SoBRLDatabaseSource::STALE_SOURCE;
    if (field == &source.inputsRevision) {
	change.reason = source.realizedInputsRevision.getValue() != source.inputsRevision.getValue() ?
	    SoBRLDatabaseSource::STALE_INPUTS : SoBRLDatabaseSource::STALE_NONE;
    } else if (field == &source.sourceRevision) {
	change.reason = source.realizedSourceRevision.getValue() != source.sourceRevision.getValue() ?
	    SoBRLDatabaseSource::STALE_SOURCE : SoBRLDatabaseSource::STALE_NONE;
    } else if (field == &source.viewRevision) {
	change.reason = source.realizedViewRevision.getValue() != source.viewRevision.getValue() ?
	    SoBRLDatabaseSource::STALE_VIEW : SoBRLDatabaseSource::STALE_NONE;
    } else if (field == &source.lodBotThreshold) {
	change.reason = SoBRLDatabaseSource::STALE_VIEW;
    } else if (field == &source.drawMode) {
	change.reason = source_external_realization_supports_draw_mode(
	    source, source.drawMode.getValue()) ?
	    SoBRLDatabaseSource::STALE_NONE : SoBRLDatabaseSource::STALE_DRAW;
    } else if (field == &source.representationMode) {
	change.reason = SoBRLDatabaseSource::STALE_DRAW;
    } else if (field == &source.tessellationAbsTol ||
	    field == &source.tessellationRelTol || field == &source.tessellationNormTol) {
	change.reason = SoBRLDatabaseSource::STALE_TESSELLATION;
    }
    return change;
}

void
SoBRLDatabaseSource::notify(SoNotList *list)
{
    SoField *field = list ? list->getLastField() : nullptr;
    const ObservedFieldChange change = this->d->observeFields ?
	observed_field_change(*this, field) : ObservedFieldChange();
    if (change.watched && !PreparedFieldNotificationScope::consume(field)) {
	if (change.reason)
	    this->failSafeObservedFieldChange(change.reason);
	this->publishObservedFieldChange(change.reason);
    }
    inherited::notify(list);
}

void
SoBRLDatabaseSource::attachObservedFieldAuditors(void)
{
    /* A quiet container suppresses its own inherited propagation, but Obol
     * skips the field notification path entirely when a field has no direct
     * auditors.  These no-op auditors keep that path live; notify() performs
     * the synchronous publication before any auditor is scheduled. */
    const std::array<SoField *, Impl::ObservedFieldCount> fields{{
	&this->path, &this->instanceKey, &this->representationKey,
	&this->representationMode, &this->drawMode, &this->tessellationAbsTol,
	&this->tessellationRelTol, &this->tessellationNormTol,
	&this->lodBotThreshold, &this->sourceRevision, &this->inputsRevision,
	&this->viewRevision
    }};
    for (size_t i = 0; i < fields.size(); ++i) {
	auto *auditor = new SoFieldSensor;
	this->d->observedFieldAuditors[i] = auditor;
	auditor->setPriority(0);
	auditor->attach(fields[i]);
    }
}

void
SoBRLDatabaseSource::detachObservedFieldAuditors(void) noexcept
{
    for (SoFieldSensor *&auditor : this->d->observedFieldAuditors) {
	if (auditor)
	    auditor->detach();
	delete auditor;
	auditor = nullptr;
    }
}

void
SoBRLDatabaseSource::clearCompiledAssembly(void)
{
    if (this->d->compiledAssembly) {
	const int childIndex = this->findChild(this->d->compiledAssembly);
	if (childIndex >= 0)
	    this->removeChild(childIndex);
	this->d->compiledAssembly->unref();
	this->d->compiledAssembly = NULL;
    }
    this->d->compiledAssemblyDirty = TRUE;
    this->d->compiledAssemblyActive = FALSE;
    this->d->compiledAssemblyNodeId = 0;
    this->d->clearCompiledCompactEvidence();
    this->markCadBatchDirty();
}

void
SoBRLDatabaseSource::ensureCompiledAssemblyChild(void)
{
    if (this->d->compiledAssembly && this->findChild(this->d->compiledAssembly) < 0)
	this->addChild(this->d->compiledAssembly);
}

void
SoBRLDatabaseSource::markCompiledAssemblyDirty(void)
{
    this->d->compiledAssemblyDirty = TRUE;
    this->d->compiledAssemblyActive = FALSE;
    this->d->compiledAssemblyNodeId = 0;
}

void
SoBRLDatabaseSource::markCadBatchDirty(void)
{
    bobol_identity_advance(this->d->cadBatchRevision);
    this->d->cadBatchDeltas.clear();
    this->d->cadBatchDeltaEntryCount = 0;
    this->d->cadBatchDeltaFloorRevision =
	this->d->cadBatchRevision;
}

/* Multiple views acknowledge independently. A journal allocation failure
 * must preserve the already-published mutation: discard sparse history and
 * force consumers below the new floor to read authoritative state once. */
template <typename Delta>
static void
database_source_record_compact_delta(uint64_t &revision,
    uint64_t &floorRevision, size_t &retainedEntryCount,
    std::deque<Delta> &deltas, const std::vector<size_t> &entryIndices,
    Delta delta)
{
    if (entryIndices.empty())
	return;
    bobol_identity_advance(revision);
    try {
	delta.revision = revision;
	delta.entryIndices = entryIndices;
	std::sort(delta.entryIndices.begin(), delta.entryIndices.end());
	delta.entryIndices.erase(std::unique(delta.entryIndices.begin(),
	    delta.entryIndices.end()), delta.entryIndices.end());
	deltas.push_back(std::move(delta));
    } catch (const std::bad_alloc &) {
	deltas.clear();
	retainedEntryCount = 0;
	floorRevision = revision;
	return;
    }
    retainedEntryCount += deltas.back().entryIndices.size();
    static const size_t maxDeltaBatches = 256;
    static const size_t maxDeltaEntries = 65536;
    while (deltas.size() > maxDeltaBatches || retainedEntryCount > maxDeltaEntries) {
	const Delta &front = deltas.front();
	retainedEntryCount -= front.entryIndices.size();
	floorRevision = std::max(floorRevision, front.revision);
	deltas.pop_front();
    }
}

void
SoBRLDatabaseSource::markCadBatchDirty(
    const std::vector<size_t> &entryIndices)
{
    database_source_record_compact_delta(this->d->cadBatchRevision,
	this->d->cadBatchDeltaFloorRevision, this->d->cadBatchDeltaEntryCount,
	this->d->cadBatchDeltas, entryIndices, Impl::CadBatchDelta());
}

SbBool
SoBRLDatabaseSource::getCadBatchChangedEntries(
    uint64_t revision, std::vector<size_t> &entryIndices) const
{
    entryIndices.clear();
    if (revision == this->d->cadBatchRevision)
	return TRUE;
    if (!revision || revision > this->d->cadBatchRevision ||
	revision < this->d->cadBatchDeltaFloorRevision)
	return FALSE;

    for (const Impl::CadBatchDelta &delta : this->d->cadBatchDeltas) {
	if (delta.revision <= revision)
	    continue;
	entryIndices.insert(entryIndices.end(), delta.entryIndices.begin(),
	    delta.entryIndices.end());
    }
    if (entryIndices.empty())
	return FALSE;
    std::sort(entryIndices.begin(), entryIndices.end());
    entryIndices.erase(std::unique(entryIndices.begin(),
	entryIndices.end()), entryIndices.end());
    return TRUE;
}

void
SoBRLDatabaseSource::markDisplayMeshLodDirty(void)
{
    bobol_identity_advance(this->d->displayMeshLodRevision);
    this->d->displayMeshLodDeltas.clear();
    this->d->displayMeshLodDeltaEntryCount = 0;
    this->d->displayMeshLodDeltaFloorRevision =
	this->d->displayMeshLodRevision;
}

void
SoBRLDatabaseSource::markDisplayMeshLodDirty(
    const std::vector<size_t> &entryIndices, SbBool coverageInvalidated)
{
    Impl::DisplayMeshLodDelta delta;
    delta.coverageInvalidated = coverageInvalidated;
    database_source_record_compact_delta(
	this->d->displayMeshLodRevision,
	this->d->displayMeshLodDeltaFloorRevision,
	this->d->displayMeshLodDeltaEntryCount,
	this->d->displayMeshLodDeltas, entryIndices, std::move(delta));
}

void
SoBRLDatabaseSource::markDisplayMeshLodVisibilityDirty(
    const std::vector<size_t> &entryIndices)
{
    database_source_record_compact_delta(
	this->d->displayMeshLodVisibilityRevision,
	this->d->displayMeshLodVisibilityDeltaFloorRevision,
	this->d->displayMeshLodVisibilityDeltaEntryCount,
	this->d->displayMeshLodVisibilityDeltas, entryIndices,
	Impl::DisplayMeshLodDelta());
}

uint64_t
SoBRLDatabaseSource::cadBatchRevisionGet(void) const
{
    return this->d->cadBatchRevision;
}

void
SoBRLDatabaseSource::clearCompactInstanceIndex(void)
{
    if (getenv("BOBOL_LOD_TRACE_SOURCE_CONTRACT") &&
	this->d->compactIndex &&
	this->d->compactIndex->sourceMeshRequestCount > 0)
	bu_log("BObol LoD source contract clearing compact index path=%s "
	       "entries=%zu requests=%zu\n",
	       this->path.getValue().getString(),
	       this->d->compactIndex->entries.size(),
	       this->d->compactIndex->sourceMeshRequestCount);
    if (this->d->compactIndex) {
	delete this->d->previousCompactIndex;
	this->d->previousCompactIndex = this->d->compactIndex;
    }
    this->d->compactIndex = NULL;
    bobol_identity_advance(this->d->compactPopulationEpoch);
    this->revokeCompactStreamAuthority();
    this->d->compactOverviewState =
	BObolCompactOccurrenceRegistryState::OverviewState::Absent;
    this->d->compactIndexActive = FALSE;
    this->d->compactOccurrenceRegistry = FALSE;
    this->d->displayMeshLodContractRevisionValid = FALSE;
    this->markCompiledAssemblyDirty();
    this->markCadBatchDirty();
    this->markDisplayMeshLodDirty();
}

void
SoBRLDatabaseSource::discardCompactInstanceHistory(void)
{
    delete this->d->compactIndex;
    this->d->compactIndex = NULL;
    bobol_identity_advance(this->d->compactPopulationEpoch);
    this->revokeCompactStreamAuthority();
    this->d->compactOverviewState =
	BObolCompactOccurrenceRegistryState::OverviewState::Absent;
    delete this->d->previousCompactIndex;
    this->d->previousCompactIndex = NULL;
    this->d->compactIndexActive = FALSE;
    this->d->compactOccurrenceRegistry = FALSE;
    this->d->displayMeshLodContractRevisionValid = FALSE;
    this->markCompiledAssemblyDirty();
    this->markCadBatchDirty();
    this->markDisplayMeshLodDirty();
}

void
SoBRLDatabaseSource::installCompactInstanceIndex(
    std::unique_ptr<BObolCompactInstanceIndex> candidate,
    SbBool occurrenceRegistry, uint64_t handleSourceId)
{
    if (!candidate)
	return;
    const BObolCompactInstanceIndex *previous = this->d->compactIndex ?
	this->d->compactIndex : this->d->previousCompactIndex;
    BObolCompactInstanceIndex *index = candidate.get();
    compact_prepare_installation(previous, *index, *this->d);
    SbBox3f bounds;
    (void)compact_index_source_bounds(*index, bounds);

    /* Every allocation above belongs to the candidate. From here, publish
     * the complete registry and release history before notifying observers. */
    this->clearCompactInstanceIndex();
    this->d->compactIndex = candidate.release();
    this->d->compactHandleSourceId = handleSourceId;
    this->d->compactIndexActive = TRUE;
    this->d->compactOccurrenceRegistry = occurrenceRegistry;
    this->d->compactExpectedInstanceCount = index->entries.size();
    for (const BObolCompactInstanceEntry &entry : index->entries) {
	if (!entry.visible || !BU_STR_EQUAL(
		entry.shapeSummary.recordRole.getString(), "lod-overview"))
	    continue;
	this->d->compactOverviewState =
	    BObolCompactOccurrenceRegistryState::OverviewState::Visible;
	break;
    }
    this->d->displayMeshLodContractRevisionValid =
	index->sourceMeshRequestCount > 0 ? TRUE : FALSE;
    this->d->displayMeshLodContractSourceRevision =
	this->sourceRevision.getValue();
    this->d->displayMeshLodContractInputsRevision =
	this->inputsRevision.getValue();
    if (getenv("BOBOL_LOD_TRACE_SOURCE_CONTRACT"))
	bu_log("BObol LoD source contract installed compact index path=%s "
	       "entries=%zu requests=%zu resident_progressive=%zu "
	       "targets=%zu registry=%d draw=%d threshold=%u "
	       "view_dependent=%d mesh_lod=%d\n",
	       this->path.getValue().getString(), index->entries.size(),
	       index->sourceMeshRequestCount,
	       index->residentProgressiveGeometryCount,
	       index->displayLodTargetCount,
	       occurrenceRegistry ? 1 : 0,
	       source_record_draw_mode(this),
	       this->lodBotThreshold.getValue(),
	       this->realizationViewDependent.getValue() ? 1 : 0,
	       this->realizationMeshLodEnabled.getValue() ? 1 : 0);
    delete this->d->previousCompactIndex;
    this->d->previousCompactIndex = NULL;

    const bool clearSelection = !index->entries.empty() &&
	!this->d->compactSelectedPaths.empty() && this->selected.getValue();
    const SbBool selectionNotifications = this->selected.enableNotify(FALSE);
    if (clearSelection)
	this->selected = FALSE;
    this->selected.enableNotify(selectionNotifications);
    if (!this->hasExactSourceBounds()) {
	if (!bounds.isEmpty())
	    (void)this->setSourceBoundsState(TRUE, bounds.getMin(), bounds.getMax(), FALSE);
	else
	    this->clearSourceBounds();
    }
    if (clearSelection && selectionNotifications)
	this->selected.touch();
}

void
SoBRLDatabaseSource::markStale(void)
{
    this->markStale(STALE_SOURCE);
}


int
SoBRLDatabaseSource::getEffectiveLodDrawMode(void) const
{
    switch (this->representationMode.getValue()) {
	case REPRESENTATION_SHADED_BOTS:
	    return BOBOL_LOD_DRAW_SHADED_BOTS;
	case REPRESENTATION_SHADED:
	    return BOBOL_LOD_DRAW_SHADED;
	case REPRESENTATION_HIDDEN_LINE:
	    return BOBOL_LOD_DRAW_HIDDEN_LINE;
	case REPRESENTATION_EVAL_POINTS:
	    return BOBOL_LOD_DRAW_POINTS;
	case REPRESENTATION_EVAL_WIRE:
	case REPRESENTATION_WIRE:
	    return BOBOL_LOD_DRAW_WIRE;
	default:
	    break;
    }

    return this->drawMode.getValue() == SHADED ?
	BOBOL_LOD_DRAW_SHADED : BOBOL_LOD_DRAW_WIRE;
}

struct SourceRealizationValues {
    std::array<uint32_t, 4> revisions;
    int status;
    /* Borrowed only during preparation; every changed field owns its copy
     * before any write can invalidate an aliased diagnostic argument. */
    const char *diagnostic;
    SbBool stale;
    uint32_t staleReason;
};

/* Only the realization record is propagated. Placement, material, compact
 * memberships and nested source owners belong to their own setters. */
template <typename StatusField, size_t AdditionalFields = 0>
class PreparedRealizationFields {
    static constexpr size_t recordFieldCount = 8;
public:
    PreparedRealizationFields(SoNode &node,
	std::array<SoSFUInt32 *, 4> revisionFields, StatusField &statusField,
	SoSFString &diagnosticField, SoSFBool &staleField,
	SoSFUInt32 &reasonField, const SourceRealizationValues &values,
	const std::array<PublicationFieldChange, AdditionalFields> &additional = {}) :
	revisions(revisionFields), status(statusField), diagnostic(diagnosticField),
	stale(staleField), reason(reasonField), next(values),
	diagnosticChanged(!database_source_string_equal(diagnosticField.getValue(), values.diagnostic)),
	preparedDiagnostic(this->diagnosticChanged ? values.diagnostic : ""),
	notifications(node, [&] {
	    std::array<PublicationFieldChange, recordFieldCount + AdditionalFields> fields{{{revisionFields[0], revisionFields[0]->getValue() != values.revisions[0]},
	    {revisionFields[1], revisionFields[1]->getValue() != values.revisions[1]},
	    {revisionFields[2], revisionFields[2]->getValue() != values.revisions[2]},
	    {revisionFields[3], revisionFields[3]->getValue() != values.revisions[3]},
	    {&statusField, statusField.getValue() != values.status},
	    {&diagnosticField, this->diagnosticChanged},
	    {&staleField, staleField.getValue() != values.stale},
	    {&reasonField, reasonField.getValue() != values.staleReason}}};
	    std::copy(additional.begin(), additional.end(), fields.begin() + recordFieldCount);
	    return fields;
	}())
    {}

    void commit()
    {
	for (size_t i = 0; i < this->revisions.size(); ++i)
	    *this->revisions[i] = this->next.revisions[i];
	this->status = this->next.status;
	if (this->diagnosticChanged)
	    this->diagnostic.setValue(std::move(this->preparedDiagnostic));
	this->stale = this->next.stale;
	this->reason = this->next.staleReason;
    }
    bool changed() const { return this->notifications.changed(); }
    void restore() { this->notifications.restore(); }
    void notify(std::exception_ptr &failure) { this->notifications.notify(failure); }

private:
    std::array<SoSFUInt32 *, 4> revisions;
    StatusField &status;
    SoSFString &diagnostic;
    SoSFBool &stale;
    SoSFUInt32 &reason;
    const SourceRealizationValues &next;
    bool diagnosticChanged;
    SbString preparedDiagnostic;
    PreparedFieldNotifications<recordFieldCount + AdditionalFields> notifications;
};

class PreparedShapeRealization {
public:
    template <typename Shape>
    PreparedShapeRealization(Shape &shape, const SourceRealizationValues &next) :
	node(shape), fields(shape, {{&shape.ownerRealizedRevision,
	    &shape.ownerRealizedSourceRevision, &shape.ownerRealizedInputsRevision,
	    &shape.ownerRealizedViewRevision}}, shape.ownerRealizationStatus,
	    shape.ownerRealizationDiagnostic, shape.ownerSourceStale,
	    shape.ownerStaleReason, next)
    {
	this->node.ref();
    }
    ~PreparedShapeRealization()
    {
	this->fields.restore();
	this->node.unref();
    }
    PreparedShapeRealization(const PreparedShapeRealization &) = delete;
    PreparedShapeRealization &operator=(const PreparedShapeRealization &) = delete;
    SoNode &node;
    PreparedRealizationFields<SoSFInt32> fields;
};

static void
prepare_shape_realization(SoNode *node, const SourceRealizationValues &next,
    std::vector<std::unique_ptr<PreparedShapeRealization>> &shapes,
    std::unordered_set<SoNode *> &visited)
{
    if (!node || node->isOfType(SoBRLDatabaseSource::getClassTypeId()) ||
	!visited.insert(node).second)
	return;
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	shapes.push_back(std::make_unique<PreparedShapeRealization>(
	    *static_cast<SoBRLVListShape *>(node), next));
    else if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	shapes.push_back(std::make_unique<PreparedShapeRealization>(
	    *static_cast<SoBRLMeshShape *>(node), next));
    else if (node->isOfType(SoGroup::getClassTypeId())) {
	auto *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); ++i)
	    prepare_shape_realization(group->getChild(i), next, shapes, visited);
    }
}

static std::vector<std::unique_ptr<PreparedShapeRealization>>
prepare_owned_realizations(SoBRLDatabaseSource &source, const SourceRealizationValues &next)
{
    std::vector<std::unique_ptr<PreparedShapeRealization>> shapes;
    std::unordered_set<SoNode *> visited;
    for (int i = 0; i < source.getNumChildren(); ++i)
	prepare_shape_realization(source.getChild(i), next, shapes, visited);
    return shapes;
}

int
SoBRLDatabaseSource::setRealizationState(int nextStatus,
	uint32_t nextRealizedSourceRevision,
	uint32_t nextRealizedInputsRevision,
	uint32_t nextStaleReason,
	const char *diagnostic)
{
    return this->publishRealizationState(nextStatus, nextRealizedSourceRevision,
	nextRealizedInputsRevision, nextStaleReason, diagnostic, nullptr);
}

int
SoBRLDatabaseSource::publishRealizationState(int nextStatus,
	uint32_t nextRealizedSourceRevision, uint32_t nextRealizedInputsRevision,
	uint32_t nextStaleReason, const char *diagnostic,
	BObolSourceRealizationEffects *effects, PublicationCommit committed, void *context)
{
    return this->publishRealizationStateWithRoles(nextStatus,
	nextRealizedSourceRevision, nextRealizedInputsRevision,
	nextStaleReason, diagnostic, this->realizationRoleFlags.getValue(),
	effects, committed, context);
}

int
SoBRLDatabaseSource::publishRealizationStateWithRoles(int nextStatus,
	uint32_t nextRealizedSourceRevision, uint32_t nextRealizedInputsRevision,
	uint32_t nextStaleReason, const char *diagnostic, int nextRoleFlags,
	BObolSourceRealizationEffects *effects, PublicationCommit committed,
	void *context)
{
    const SbBool wasBatchEligible =
	this->realizationStatus.getValue() == REALIZED &&
	!this->stale.getValue();
    if (nextStatus != REALIZED && nextStatus != UNREALIZED &&
	nextStatus != FAILED)
	nextStatus = UNREALIZED;

    const SbBool realized = nextStatus == REALIZED ? TRUE : FALSE;
    if (realized) {
	nextRealizedSourceRevision = nextRealizedSourceRevision ?
				     nextRealizedSourceRevision : this->sourceRevision.getValue();
	nextRealizedInputsRevision = nextRealizedInputsRevision ?
				     nextRealizedInputsRevision : this->inputsRevision.getValue();
	nextStaleReason = STALE_NONE;
    } else if (!nextStaleReason) {
	nextStaleReason = STALE_SOURCE;
    }

    const SbBool nextStale = realized ? FALSE : TRUE;
    const SourceRealizationValues next{{{
	realized ? nextRealizedSourceRevision : this->realizedRevision.getValue(),
	realized ? nextRealizedSourceRevision : this->realizedSourceRevision.getValue(),
	realized ? nextRealizedInputsRevision : this->realizedInputsRevision.getValue(),
	realized ? this->viewRevision.getValue() : this->realizedViewRevision.getValue()}},
	nextStatus, diagnostic ? diagnostic : "", nextStale, nextStaleReason};
    const int validRoleFlags = REALIZATION_ROLE_CSG |
	REALIZATION_ROLE_MESH | REALIZATION_ROLE_EXTERNAL;
    nextRoleFlags &= validRoleFlags;
    const bool rolesChanged =
	this->realizationRoleFlags.getValue() != nextRoleFlags;
    PreparedRealizationFields<SoSFEnum, 1> publication(*this,
	{{&this->realizedRevision, &this->realizedSourceRevision,
	    &this->realizedInputsRevision, &this->realizedViewRevision}},
	this->realizationStatus, this->realizationDiagnostic, this->stale,
	this->staleReason, next,
	{{{&this->realizationRoleFlags, rolesChanged}}});
    if (!publication.changed())
	return 0;
    auto shapes = prepare_owned_realizations(*this, next);

    publication.commit();
    if (rolesChanged)
	this->realizationRoleFlags = nextRoleFlags;
    for (auto &shape : shapes)
	shape->fields.commit();
    const SbBool isBatchEligible = realized && !nextStale;
    if (wasBatchEligible != isBatchEligible || rolesChanged)
	this->markCadBatchDirty();
    this->markCompiledAssemblyDirty();
    if (effects) effects->commit(true);
    if (committed) committed(context);
    publication.restore();
    for (auto &shape : shapes)
	shape->fields.restore();
    std::exception_ptr failure;
    publication.notify(failure);
    for (auto &shape : shapes)
	shape->fields.notify(failure);
    if (effects) effects->notify(failure);
    if (failure)
	std::rethrow_exception(failure);
    return 1;
}

static bool
invalidation_revokes_source_contract(uint32_t reason)
{
    return (reason & (SoBRLDatabaseSource::STALE_SOURCE | SoBRLDatabaseSource::STALE_INPUTS |
	SoBRLDatabaseSource::STALE_DATABASE)) != 0;
}

void
SoBRLDatabaseSource::markStale(uint32_t reason)
{
    if (!reason) return;
    const SbBool hadCurrentDisplayMeshLodContract = this->hasDisplayMeshLodRequests();
    /* Only source identity changes revoke immutable producer data. View and
     * presentation changes retain that data for the existing LoD planner. */
    const bool revokeSourceContract = invalidation_revokes_source_contract(reason);
    const bool retainExternalFailure = reason == STALE_VIEW &&
	(this->realizationRoleFlags.getValue() & REALIZATION_ROLE_EXTERNAL) &&
	this->realizationStatus.getValue() == FAILED;
    const SourceRealizationValues next{{{this->realizedRevision.getValue(),
	this->realizedSourceRevision.getValue(), this->realizedInputsRevision.getValue(),
	this->realizedViewRevision.getValue()}}, retainExternalFailure ? FAILED : UNREALIZED,
	retainExternalFailure ? this->realizationDiagnostic.getValue().getString() : "",
	TRUE, this->staleReason.getValue() | reason};
    PreparedRealizationFields<SoSFEnum, 1> publication(*this,
	{{&this->realizedRevision, &this->realizedSourceRevision,
	    &this->realizedInputsRevision, &this->realizedViewRevision}},
	this->realizationStatus, this->realizationDiagnostic, this->stale, this->staleReason,
	next, {{{&this->sourceBoundsExact, revokeSourceContract && this->sourceBoundsExact.getValue()}}});
    auto shapes = prepare_owned_realizations(*this, next);

    publication.commit();
    for (auto &shape : shapes) shape->fields.commit();
    if (revokeSourceContract) {
	this->revokeSourceResourceAuthority();
	this->sourceBoundsExact = FALSE;
    }
    this->markCompiledAssemblyDirty();
    this->markCadBatchDirty();
    publication.restore();
    for (auto &shape : shapes) shape->fields.restore();
    std::exception_ptr failure;
    publication.notify(failure);
    for (auto &shape : shapes) shape->fields.notify(failure);
    try {
	if (getenv("BOBOL_LOD_TRACE_SOURCE_CONTRACT") && hadCurrentDisplayMeshLodContract &&
	    !this->hasDisplayMeshLodRequests())
	    bu_log("BObol LoD source contract invalidated path=%s reason=%u "
		"source_revision=%u inputs_revision=%u requests=%zu\n",
		this->path.getValue().getString(), reason, this->sourceRevision.getValue(),
		this->inputsRevision.getValue(), this->d->compactIndex ?
		    this->d->compactIndex->sourceMeshRequestCount : 0);
    } catch (...) { if (!failure) failure = std::current_exception(); }
    if (failure) std::rethrow_exception(failure);
}

int
SoBRLDatabaseSource::setRealizationRoleFlags(int roleFlags)
{
    return this->setRealizationRoleFlags(roleFlags, nullptr, nullptr);
}

int
SoBRLDatabaseSource::setRealizationRoleFlags(int roleFlags, PublicationCommit committed, void *context)
{
    const int validFlags = REALIZATION_ROLE_CSG | REALIZATION_ROLE_MESH |
			   REALIZATION_ROLE_EXTERNAL;
    roleFlags &= validFlags;
    if (this->realizationRoleFlags.getValue() == roleFlags)
	return 0;

    PreparedFieldNotifications<1> publication(*this, {{{&this->realizationRoleFlags, true}}});
    this->realizationRoleFlags = roleFlags;
    // Role flags affect source admission; owned geometry and appearance stay intact.
    this->markCompiledAssemblyDirty();
    this->markCadBatchDirty();
    if (committed) committed(context);
    std::exception_ptr failure;
    publication.notify(failure);
    if (failure) std::rethrow_exception(failure);
    return 1;
}

int
SoBRLDatabaseSource::setSourceBoundsState(SbBool nextBoundsValid,
	const SbVec3f &nextBoundsMin,
	const SbVec3f &nextBoundsMax,
	SbBool nextBoundsExact)
{
    return this->setSourceBoundsState(nextBoundsValid, nextBoundsMin, nextBoundsMax,
	nextBoundsExact, nullptr, nullptr);
}

int
SoBRLDatabaseSource::setSourceBoundsState(SbBool nextBoundsValid,
	const SbVec3f &nextBoundsMin, const SbVec3f &nextBoundsMax,
	SbBool nextBoundsExact, PublicationCommit committed, void *context)
{
    SbVec3f sanitizedMin(0.0f, 0.0f, 0.0f);
    SbVec3f sanitizedMax(0.0f, 0.0f, 0.0f);
    if (nextBoundsValid) {
	const SbBox3f bounds =
	    database_source_box_from_minmax(nextBoundsMin, nextBoundsMax);
	sanitizedMin = bounds.getMin();
	sanitizedMax = bounds.getMax();
    }

    nextBoundsExact = nextBoundsValid && nextBoundsExact;
    const std::array<SoField *, 4> fields = {{&this->sourceBoundsValid,
	&this->sourceBoundsExact, &this->sourceBoundsMin, &this->sourceBoundsMax}};
    const std::array<bool, fields.size()> differences = {{
	this->sourceBoundsValid.getValue() != nextBoundsValid,
	this->sourceBoundsExact.getValue() != nextBoundsExact,
	!database_source_vec3f_equal(this->sourceBoundsMin.getValue(), sanitizedMin),
	!database_source_vec3f_equal(this->sourceBoundsMax.getValue(), sanitizedMax)}};
    std::array<SbBool, fields.size()> notifications;
    for (size_t i = 0; i < fields.size(); ++i)
	notifications[i] = fields[i]->enableNotify(FALSE);
    /* Observers must see one coherent bound even if notification allocation
     * fails. Numeric field assignment with notifications disabled allocates
     * nothing; restore every notification flag before invoking any observer. */
    if (differences[0])
	this->sourceBoundsValid = nextBoundsValid;
    if (differences[1])
	this->sourceBoundsExact = nextBoundsExact;
    if (differences[2])
	this->sourceBoundsMin = sanitizedMin;
    if (differences[3])
	this->sourceBoundsMax = sanitizedMax;
    for (size_t i = 0; i < fields.size(); ++i)
	fields[i]->enableNotify(notifications[i]);
    const bool changed = std::any_of(differences.begin(), differences.end(), [](bool value) { return value; });
    if (changed && committed) committed(context);
    std::exception_ptr failure;
    for (size_t i = 0; i < fields.size(); ++i) {
	if (!differences[i] || !notifications[i])
	    continue;
	try { fields[i]->touch(); }
	catch (...) { if (!failure) failure = std::current_exception(); }
    }
    if (failure) std::rethrow_exception(failure);
    return changed;
}

int
SoBRLDatabaseSource::setSourceBoundsExactState(SbBool nextBoundsExact)
{
    nextBoundsExact = nextBoundsExact && this->sourceBoundsValid.getValue();
    if (this->sourceBoundsExact.getValue() == nextBoundsExact)
	return 0;
    this->sourceBoundsExact = nextBoundsExact;
    return 1;
}

void
SoBRLDatabaseSource::clearSourceBounds(void)
{
    (void)this->setSourceBoundsState(FALSE,
				     SbVec3f(0.0f, 0.0f, 0.0f),
				     SbVec3f(0.0f, 0.0f, 0.0f), FALSE);
}

SbBool
SoBRLDatabaseSource::getSourceBounds(SbBox3f &bounds) const
{
    bounds.makeEmpty();
    if (!this->sourceBoundsValid.getValue())
	return FALSE;

    bounds = database_source_box_from_minmax(
		 this->sourceBoundsMin.getValue(),
		 this->sourceBoundsMax.getValue());
    return bounds.isEmpty() ? FALSE : TRUE;
}

SbBool
SoBRLDatabaseSource::getEffectiveSourceBounds(SbBox3f &bounds) const
{
    if (!this->getSourceBounds(bounds))
	return FALSE;

    const SoMatrixTransform *transform = source_placement_transform(this);
    if (transform) {
	bounds = database_source_transform_bounds(bounds,
		 transform->matrix.getValue());
    } else if (this->drawMatrixValid.getValue()) {
	bounds = database_source_transform_bounds(bounds,
		 this->drawMatrix.getValue());
    }

    return bounds.isEmpty() ? FALSE : TRUE;
}

SbBool
SoBRLDatabaseSource::hasExactSourceBounds(void) const
{
    return this->sourceBoundsValid.getValue() &&
	this->sourceBoundsExact.getValue();
}

bool
BObolCadSourceState::setDatabaseBinding(struct db_i *database)
{
    if (this->dbip == database)
	return false;
    this->dbip = database;
    bobol_identity_advance(this->databaseBindingEpoch);
    return true;
}

void
SoBRLDatabaseSource::setDatabase(struct db_i *database)
{
    if (this->d->setDatabaseBinding(database))
	this->markStale(STALE_DATABASE);
}

struct db_i *
SoBRLDatabaseSource::getDatabase(void) const {
    return this->d->dbip;
}

static bool
database_source_mesh_lod_bounds_valid(const SbVec3f &bmin,
	const SbVec3f &bmax)
{
    for (int axis = 0; axis < 3; ++axis) {
	if (!std::isfinite(bmin[axis]) || !std::isfinite(bmax[axis]) ||
		bmin[axis] > bmax[axis])
	    return false;
    }
    return true;
}

int
SoBRLDatabaseSource::adoptMeshLod(
    const BObolSourceRealizationStamp &stamp,
    struct BObolMeshLod *lod,
    const SbVec3f &bmin,
    const SbVec3f &bmax)
{
    if (!lod || lod == this->d->meshLod ||
	!this->matchesRealizationStamp(stamp) ||
	!database_source_mesh_lod_bounds_valid(bmin, bmax))
	return 0;

    struct BObolMeshLod *retired = this->d->meshLod;
    this->d->meshLod = lod;
    this->d->meshLodBoundsMin = bmin;
    this->d->meshLodBoundsMax = bmax;
    this->d->meshLodBoundsValid = TRUE;
    if (retired)
	bobol_mesh_lod_destroy(retired);
    return 1;
}

struct BObolMeshLod *
SoBRLDatabaseSource::getMeshLod(void) const {
    return this->d->meshLod;
}

void
SoBRLDatabaseSource::clearMeshLod(void)
{
    struct BObolMeshLod *retired = this->d->meshLod;
    this->d->meshLod = NULL;
    this->d->meshLodBoundsValid = FALSE;
    this->d->meshLodBoundsMin.setValue(0.0f, 0.0f, 0.0f);
    this->d->meshLodBoundsMax.setValue(0.0f, 0.0f, 0.0f);
    if (retired)
	bobol_mesh_lod_destroy(retired);
}

SbBool
SoBRLDatabaseSource::getMeshLodBounds(SbVec3f &bmin,
				      SbVec3f &bmax) const
{
    if (!this->d->meshLodBoundsValid)
	return FALSE;

    bmin = this->d->meshLodBoundsMin;
    bmax = this->d->meshLodBoundsMax;
    return TRUE;
}

void
SoBRLDatabaseSource::revokeCompactStreamAuthority(void)
{
    this->d->compactExpectedInstanceCount = 0;
    this->d->compactExpectedInstanceCountCertified = FALSE;
    this->d->compactSourceProfile = BObolCompactSourceProfile();
    this->d->compactStagedSourceStream.reset();
}

void
SoBRLDatabaseSource::revokeSourceResourceAuthority(void)
{
    this->d->displayMeshLodContractRevisionValid = FALSE;
    this->revokeCompactStreamAuthority();
    this->clearMeshLod();
}

void
SoBRLDatabaseSource::configureDatabaseSource(const char *sourcePath,
	struct db_i *database,
	int mode,
	uint32_t revision)
{
    this->configureDatabaseSourceInstance(sourcePath, sourcePath, database,
					  mode, revision);
}

void
SoBRLDatabaseSource::configureDatabaseSourceInstance(
    const char *sourceInstanceKey,
    const char *sourcePath,
    struct db_i *database,
    int mode,
    uint32_t revision)
{
    this->configureDatabaseSourceInstanceRepresentation(sourceInstanceKey,
	    sourcePath, NULL, -1, database, mode, revision);
}

static void
database_snapshot_collect_tree_dependencies(const union tree *tree,
	std::vector<std::string> &dependencies)
{
    if (!tree)
	return;

    switch (tree->tr_op) {
	case OP_UNION:
	case OP_INTERSECT:
	case OP_SUBTRACT:
	case OP_XOR:
	    database_snapshot_collect_tree_dependencies(tree->tr_b.tb_right,
		dependencies);
	    /* fall through */
	case OP_NOT:
	case OP_GUARD:
	case OP_XNOP:
	    database_snapshot_collect_tree_dependencies(tree->tr_b.tb_left,
		dependencies);
	    break;
	case OP_DB_LEAF:
	    if (tree->tr_l.tl_name && tree->tr_l.tl_name[0])
		dependencies.push_back(tree->tr_l.tl_name);
	    break;
	default:
	    break;
    }
}

static int
database_snapshot_copy_object(struct db_i *source,
	struct rt_wdb *targetWdb, const char *name,
	std::unordered_set<std::string> &copied,
	std::unordered_set<std::string> &visiting)
{
    if (!source || !targetWdb || !name || !name[0])
	return 0;

    const std::string objectName(name);
    if (copied.find(objectName) != copied.end())
	return 1;
    if (visiting.find(objectName) != visiting.end())
	return 1;

    struct directory *dp = db_lookup(source, name, LOOKUP_QUIET);
	if (!dp) {
	bu_log("Obol detached realization snapshot could not find '%s'\n", name);
	return 0;
	}
    visiting.insert(objectName);

    if (dp->d_flags & RT_DIR_COMB) {
	struct rt_db_internal intern;
	RT_DB_INTERNAL_INIT(&intern);
	if (rt_db_get_internal(&intern, dp, source, NULL) < 0) {
	    bu_log("Obol detached realization snapshot could not read combination '%s'\n",
		name);
	    visiting.erase(objectName);
	    return 0;
	}
	std::vector<std::string> dependencies;
	if (intern.idb_type == ID_COMBINATION && intern.idb_ptr) {
	    const struct rt_comb_internal *comb =
		static_cast<const struct rt_comb_internal *>(intern.idb_ptr);
	    database_snapshot_collect_tree_dependencies(comb->tree,
		dependencies);
	}
	rt_db_free_internal(&intern);
	for (const std::string &dependency : dependencies) {
	    if (!database_snapshot_copy_object(source, targetWdb,
		    dependency.c_str(), copied, visiting)) {
		bu_log("Obol detached realization snapshot could not copy dependency '%s' of '%s'\n",
		    dependency.c_str(), name);
		visiting.erase(objectName);
		return 0;
	    }
	}
    }

    struct bu_external external;
    BU_EXTERNAL_INIT(&external);
    if (db_get_external(&external, dp, source) < 0) {
	bu_log("Obol detached realization snapshot could not read '%s'\n", name);
	visiting.erase(objectName);
	return 0;
    }
    const int copiedObject = wdb_export_external(targetWdb, &external,
	name, dp->d_flags & ~RT_DIR_INMEM, dp->d_minor_type) >= 0;
    bu_free_external(&external);
    visiting.erase(objectName);
    if (!copiedObject)
	{
	bu_log("Obol detached realization snapshot could not export '%s'\n", name);
	return 0;
	}

    copied.insert(objectName);
    return 1;
}

static std::string
database_snapshot_root_name(const SbString &sourcePath)
{
    const std::string path = database_lookup_path_from_source_path(sourcePath);
    const std::string::size_type slash = path.find('/');
    return slash == std::string::npos ? path : path.substr(0, slash);
}

static void
database_snapshot_copy_lookup_context(struct db_i *snapshot,
	const struct db_i *source)
{
    if (!snapshot || !source)
	return;

    /* Compact progressive realization identifies each combination leaf by
     * occurrence, not merely by the referenced directory object.  Detached
     * database snapshots must therefore walk with instance specifiers enabled
     * just like the live draw path.  Depending on the process environment here
     * made cold and warm draws disagree, and duplicate instances could collapse
     * into an unusable empty compact realization. */
    (void)db_comb_instance_ids_set(snapshot, 1);
    snapshot->dbi_local2base = source->dbi_local2base;
    snapshot->dbi_base2local = source->dbi_base2local;
    if (snapshot->dbi_filename) {
	bu_free(snapshot->dbi_filename, "database snapshot filename");
	snapshot->dbi_filename = NULL;
    }
    if (source->dbi_filename)
	snapshot->dbi_filename = bu_strdup(source->dbi_filename);
    if (snapshot->dbi_title) {
	bu_free(snapshot->dbi_title, "database snapshot title");
	snapshot->dbi_title = NULL;
    }
    if (source->dbi_title)
	snapshot->dbi_title = bu_strdup(source->dbi_title);
    if (snapshot->dbi_filepath) {
	bu_argv_free(2, snapshot->dbi_filepath);
	snapshot->dbi_filepath = NULL;
    }
    if (!source->dbi_filepath)
	return;

    /* db_close owns the standard two search entries, so retain that ABI even
     * for an in-memory snapshot. */
    snapshot->dbi_filepath = static_cast<char **>(bu_calloc(3,
	sizeof(char *), "database snapshot filepath"));
    snapshot->dbi_filepath[0] = bu_strdup(source->dbi_filepath[0] ?
	source->dbi_filepath[0] : ".");
    snapshot->dbi_filepath[1] = bu_strdup(source->dbi_filepath[1] ?
	source->dbi_filepath[1] : ".");
}

static struct db_i *
database_snapshot_create(struct db_i *source, const SbString &sourcePath,
	SbString *snapshotPathOut)
{
    if (snapshotPathOut)
	*snapshotPathOut = "";
    if (!source)
	return NULL;

    const std::string rootName = database_snapshot_root_name(sourcePath);
    const int version = db_version(source);
    struct db_i *snapshot = NULL;
    SbString snapshotPath;
    int wdbType = RT_WDB_TYPE_DB_INMEM;

    if (rootName.empty() || (version != 4 && version != 5))
	return NULL;
    if (version == 5) {
	snapshot = db_open_inmem();
    } else {
	char path[MAXPATHLEN] = {0};
	FILE *file = bu_temp_file(path, sizeof(path));
	if (!file)
	    return NULL;
	(void)fclose(file);
	(void)bu_file_delete(path);
	snapshot = db_create(path, version);
	if (snapshot)
	    snapshotPath = path;
	wdbType = RT_WDB_TYPE_DB_DISK;
    }
    if (!snapshot)
	{
	bu_log("Obol detached realization could not create a v%d database snapshot\n",
	    version);
	return NULL;
	}

    struct rt_wdb *snapshotWdb = wdb_dbopen(snapshot, wdbType);
    std::unordered_set<std::string> copied;
    std::unordered_set<std::string> visiting;
    int success = snapshotWdb ? 1 : 0;
    struct directory *global = db_lookup(source, DB5_GLOBAL_OBJECT_NAME,
	LOOKUP_QUIET);
    if (success && global)
	success = database_snapshot_copy_object(source, snapshotWdb,
	    DB5_GLOBAL_OBJECT_NAME, copied, visiting);
    if (success)
	success = database_snapshot_copy_object(source, snapshotWdb,
	    rootName.c_str(), copied, visiting);
    if (!success) {
	bu_log("Obol detached realization could not copy source closure '%s' into a v%d database snapshot\n",
	    rootName.c_str(), version);
	db_close(snapshot);
	if (snapshotPath.getLength() > 0)
	    (void)bu_file_delete(snapshotPath.getString());
	return NULL;
    }

    if (wdbType == RT_WDB_TYPE_DB_DISK)
	db_sync(snapshot);
    /* Copying _GLOBAL writes the persistent attribute, but an in-memory
     * snapshot is not rescanned afterward.  Populate its runtime material
     * table explicitly so detached progressive realization uses the same
     * active region-id colors as the live database. */
    struct bu_vls colorTable = BU_VLS_INIT_ZERO;
    db_mater_to_vls(&colorTable, source);
    if (bu_vls_strlen(&colorTable) > 0)
	db5_import_color_table(snapshot, bu_vls_addr(&colorTable));
    bu_vls_free(&colorTable);
    database_snapshot_copy_lookup_context(snapshot, source);
    snapshot->dbi_read_only = 1;
    if (snapshotPathOut)
	*snapshotPathOut = snapshotPath;
    return snapshot;
}

SbBool
BObolSourceRealizationStamp::matchesConfiguration(
    const SoBRLDatabaseSource &source) const
{
    return this->instanceKey == source.instanceKey.getValue() &&
	this->path == source.path.getValue() &&
	this->representationKey == source.representationKey.getValue() &&
	this->drawMode == source.drawMode.getValue() &&
	this->representationMode == source.representationMode.getValue() &&
	this->sourceRevision == source.sourceRevision.getValue() &&
	this->inputsRevision == source.inputsRevision.getValue() &&
	this->tessellationTolerance == SbVec3f(source.tessellationAbsTol.getValue(),
	    source.tessellationRelTol.getValue(), source.tessellationNormTol.getValue());
}

BObolSourceRealizationStamp
SoBRLDatabaseSource::captureRealizationStamp(void) const
{
    BObolSourceRealizationStamp stamp;
    stamp.routingId = this->d->routingId;
    stamp.populationEpoch = this->d->compactPopulationEpoch;
    stamp.databaseBindingEpoch = this->d->databaseBindingEpoch;
    stamp.instanceKey = this->instanceKey.getValue();
    stamp.path = this->path.getValue();
    stamp.representationKey = this->representationKey.getValue();
    stamp.drawMode = this->drawMode.getValue();
    stamp.representationMode = this->representationMode.getValue();
    stamp.sourceRevision = this->sourceRevision.getValue();
    stamp.inputsRevision = this->inputsRevision.getValue();
    stamp.tessellationTolerance.setValue(this->tessellationAbsTol.getValue(),
	this->tessellationRelTol.getValue(), this->tessellationNormTol.getValue());
    return stamp;
}

SbBool
SoBRLDatabaseSource::matchesRealizationStamp(
    const BObolSourceRealizationStamp &stamp) const
{
    return stamp.routingId == this->d->routingId &&
	stamp.populationEpoch == this->d->compactPopulationEpoch &&
	stamp.databaseBindingEpoch == this->d->databaseBindingEpoch &&
	stamp.matchesConfiguration(*this);
}

BObolSourcePresentationStamp
SoBRLDatabaseSource::capturePresentationStamp(void) const
{
    BObolSourcePresentationStamp stamp;
    stamp.routingId = this->d->routingId;
    stamp.cadRevision = this->cadBatchRevisionGet();
    return stamp;
}

SbBool
SoBRLDatabaseSource::matchesPresentationStamp(
    const BObolSourcePresentationStamp &stamp) const
{
    return stamp.routingId == this->d->routingId &&
	stamp.cadRevision == this->cadBatchRevisionGet();
}

SoBRLDatabaseSource *
SoBRLDatabaseSource::createDetachedRealizationTemplate(void) const
{
    SbModernUtils::SoNodeRef owner(
	new SoBRLDatabaseSource(FieldObservation::Detached));
    auto *detached = static_cast<SoBRLDatabaseSource *>(owner.get());
    detached->copyFieldValues(this, FALSE);
    detached->stale = TRUE;
    detached->staleReason = STALE_SOURCE;
    detached->realizationStatus = UNREALIZED;
    detached->realizationDiagnostic = "";
    /* Private construction uses immutable configuration, either locally or
     * exclusively on a worker. It needs no owner-thread invalidation sensors;
     * only its completed outputs become observable at publication. */
    owner.release();
    return detached;
}

SbBool
SoBRLDatabaseSource::initializeDetachedRealizationDatabase(
	struct db_i *sourceDatabase, struct db_i **databaseOut,
	SbString *snapshotPathOut)
{
    if (databaseOut)
	*databaseOut = NULL;
    if (snapshotPathOut)
	*snapshotPathOut = "";
    if (!databaseOut || !sourceDatabase || this->d->dbip)
	return FALSE;

    /*
     * A file-backed database already has an immutable directory index and a
     * librt-managed, serialized I/O path.  Retain that database instance for
     * the worker instead of reopening the file and rebuilding its complete
     * directory before the first coverage box can be published.  The latter
     * made every cold draw O(database object count), even when the requested
     * closure was small, and dominated time-to-first-pixel on vehicle-scale
     * databases.
     *
     * db_clone_dbi() is librt's ordinary additional-client contract (also
     * used by rt_i and the LoD database leases).  It retains the indexed
     * database while db_read() serializes access to shared file state.  This
     * detached source is read-only by contract; source/input revisions still
     * reject streamed or final results if an edit supersedes the request.
     *
     * An in-memory database still needs a closure snapshot: its records can
     * be replaced in place and there is no persistent file/index to retain.
     */
    struct db_i *database = NULL;
    if (sourceDatabase->dbi_filename && sourceDatabase->dbi_filename[0])
	database = db_clone_dbi(sourceDatabase, NULL);
    if (!database)
	database = database_snapshot_create(sourceDatabase,
	    this->path.getValue(), snapshotPathOut);
    if (!database)
	return FALSE;
    (void)this->d->setDatabaseBinding(database);
    *databaseOut = database;
    return TRUE;
}

SoBRLDatabaseSource *
SoBRLDatabaseSource::createDetachedRealizationSource(
	struct db_i **databaseOut, SbString *snapshotPathOut) const
{
    if (databaseOut)
	*databaseOut = NULL;
    if (snapshotPathOut)
	*snapshotPathOut = "";
    if (!databaseOut || !this->d->dbip)
	return NULL;

    SoBRLDatabaseSource *detached =
	this->createDetachedRealizationTemplate();
    if (!detached)
	return NULL;
    SbModernUtils::SoNodeRef owner(detached);
    detached->unref();
    struct db_i *database = database_snapshot_create(this->d->dbip,
	this->path.getValue(), snapshotPathOut);
    if (!database)
	return NULL;
    (void)detached->d->setDatabaseBinding(database);
    *databaseOut = database;
    owner.release();
    return detached;
}

/* Prepare strings and notification state before changing the source. The fixed
 * field set is the terminal realization record, not another scene snapshot. */
class PreparedSourceRealization {
public:
    PreparedSourceRealization(SoBRLDatabaseSource &target,
	const SbBox3f &bounds, SbBool exact, bool clearRootSelection, int roleFlags) :
	PreparedSourceRealization(target, target.sourceRevision.getValue(), bounds, exact, clearRootSelection, roleFlags)
    {}

    PreparedSourceRealization(SoBRLDatabaseSource &target, uint32_t nextRevision,
	const SbBox3f &bounds, SbBool exact, bool clearRootSelection, int roleFlags) :
	source(target), revision(nextRevision),
	boundsValid(!bounds.isEmpty()),
	boundsExact(this->boundsValid && exact),
	minimum(this->boundsValid ? bounds.getMin() : SbVec3f(0, 0, 0)),
	maximum(this->boundsValid ? bounds.getMax() : SbVec3f(0, 0, 0)),
	identity(source_realization_identity(&target, nextRevision)),
	clearSelection(clearRootSelection),
	roles(roleFlags),
	notifications(target, {{
	    {&target.sourceRevision, target.sourceRevision.getValue() != this->revision},
	    {&target.realizedRevision, target.realizedRevision.getValue() != this->revision},
	    {&target.realizedSourceRevision, target.realizedSourceRevision.getValue() != this->revision},
	    {&target.realizedInputsRevision, target.realizedInputsRevision.getValue() != target.inputsRevision.getValue()},
	    {&target.realizedViewRevision, target.realizedViewRevision.getValue() != target.viewRevision.getValue()},
	    {&target.realizationStatus, target.realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED},
	    {&target.realizationDiagnostic, target.realizationDiagnostic.getValue().getLength() != 0},
	    {&target.realizationIdentity, (target.realizationIdentity.getValue() != this->identity) != FALSE},
	    {&target.stale, target.stale.getValue() != FALSE},
	    {&target.staleReason, target.staleReason.getValue() != SoBRLDatabaseSource::STALE_NONE},
	    {&target.sourceBoundsValid, target.sourceBoundsValid.getValue() != this->boundsValid},
	    {&target.sourceBoundsExact, target.sourceBoundsExact.getValue() != this->boundsExact},
	    {&target.sourceBoundsMin, !database_source_vec3f_equal(target.sourceBoundsMin.getValue(), this->minimum)},
	    {&target.sourceBoundsMax, !database_source_vec3f_equal(target.sourceBoundsMax.getValue(), this->maximum)},
	    {&target.selected, clearRootSelection && target.selected.getValue()},
	    {&target.realizationRoleFlags, target.realizationRoleFlags.getValue() != roleFlags}
	}})
    {}

    PreparedSourceRealization(const PreparedSourceRealization &) = delete;
    PreparedSourceRealization &operator=(const PreparedSourceRealization &) = delete;

    void commit()
    {
	/* String storage transfers into a quiet field. Numeric setters and
	 * shared ownership publication likewise need no further allocation. */
	(void)this->source.setSourceBoundsState(this->boundsValid,
	    this->minimum, this->maximum, this->boundsExact);
	this->source.sourceRevision = this->revision;
	this->source.realizedRevision = this->source.sourceRevision.getValue();
	this->source.realizedSourceRevision = this->source.sourceRevision.getValue();
	this->source.realizedInputsRevision = this->source.inputsRevision.getValue();
	this->source.realizedViewRevision = this->source.viewRevision.getValue();
	this->source.realizationStatus = SoBRLDatabaseSource::REALIZED;
	this->source.realizationDiagnostic.setValue(SbString());
	this->source.realizationIdentity.setValue(std::move(this->identity));
	this->source.stale = FALSE;
	this->source.staleReason = SoBRLDatabaseSource::STALE_NONE;
	this->source.realizationRoleFlags = this->roles;
	if (this->clearSelection)
	    this->source.selected = FALSE;
    }

    void restore() { this->notifications.restore(); }
    void notify(std::exception_ptr &failure) { this->notifications.notify(failure); }
    void notify()
    {
	std::exception_ptr failure;
	this->notify(failure);
	if (failure)
	    std::rethrow_exception(failure);
    }

private:
    SoBRLDatabaseSource &source;
    uint32_t revision;
    SbBool boundsValid;
    SbBool boundsExact;
    SbVec3f minimum;
    SbVec3f maximum;
    SbString identity;
    bool clearSelection;
    int roles;
    PreparedFieldNotifications<16> notifications;
};

int
SoBRLDatabaseSource::adoptDetachedCompactRealization(
    const BObolSourceRealizationStamp &stamp,
    SoBRLDatabaseSource *detached,
    SbBool authoritativeStreamDrained,
    const std::shared_ptr<BObolCompactOccurrenceStream> &stagedSourceStream)
{
    return this->adoptDetachedCompactRealization(stamp, detached,
	authoritativeStreamDrained, stagedSourceStream, nullptr);
}

int
SoBRLDatabaseSource::adoptDetachedCompactRealization(
    const BObolSourceRealizationStamp &stamp,
    SoBRLDatabaseSource *detached,
    SbBool authoritativeStreamDrained,
    const std::shared_ptr<BObolCompactOccurrenceStream> &stagedSourceStream,
    BObolSourceRealizationEffects *effects)
{
    if (!detached || !this->matchesRealizationStamp(stamp) ||
	!stamp.matchesConfiguration(*detached))
	return 0;

    SbModernUtils::SoNodeRef sourceOwner(this);
    SbModernUtils::SoNodeRef detachedOwner(detached);

    /*
     * A progressive worker publishes every authoritative occurrence through
     * the hand-off stream before it reaches COMPLETE.  When that stream has
     * been drained, the live index already is the detached realization; its
     * only possible extra is the temporary whole-target overview.  Preserve
     * the live index and compiled assembly in that common case.  Replacing
     * the index here used to invalidate every retained part generation and
     * turn the final hand-off of a 5k/50k draw into one unbounded rebuild.
     */
    BObolCompactInstanceIndex *current = this->d->compactIndex;
    BObolCompactInstanceIndex *authoritative = detached->d->compactIndex;
    if (!current || current->entries.empty() ||
	(!authoritative && !authoritativeStreamDrained))
	return 0;
    size_t overviewCount = 0;
    for (BObolCompactInstanceEntry &entry : current->entries) {
	if (BU_STR_EQUAL(entry.shapeSummary.recordRole.getString(),
		"lod-overview"))
	    overviewCount++;
    }
    const size_t streamedAuthoritativeCount =
	current->entries.size() - overviewCount;
    const size_t authoritativeCount = authoritative ?
	authoritative->entries.size() : streamedAuthoritativeCount;
    /* An empty detached registry is not authority to remove the only useful
     * cold-start presentation.  This happens for deferred roots whose leaf
     * producer finishes its structural bookkeeping before a drawable leaf has
     * reached the stream.  Retiring the overview at that point leaves a blank
     * view until a later LoD request happens to publish. */
    /* Only the caller that merged the complete authoritative stream can
     * certify reuse.  Matching paths, tiers and transforms can still name
     * different immutable geometry or provider requests after a source edit.
     * Without that certificate, install the detached registry below. */
    const bool streamedComplete = authoritativeStreamDrained &&
	streamedAuthoritativeCount > 0 &&
	(!authoritative || streamedAuthoritativeCount == authoritativeCount);
    if (!streamedComplete && (!authoritative || authoritative->entries.empty()))
	return 0;
    if (bobol_transaction_fault_requested(
	    BObolTransactionFaultPoint::SOURCE_TERMINAL_PREPARATION))
	throw std::bad_alloc();

    SbBox3f bounds;
    SbBool boundsExact = FALSE;
    SbBox3f streamBounds;
    const bool haveStreamBounds = stagedSourceStream &&
	stagedSourceStream->getCoverageBounds(streamBounds) &&
	!streamBounds.isEmpty();
    /* Producer-certified whole-root coverage survives display-derived bounds
     * and belongs to this terminal publication, not a later setter call. */
    if (haveStreamBounds) {
	bounds = streamBounds;
	boundsExact = stagedSourceStream->hasCoverageBoundsComplete() ? TRUE : FALSE;
    } else if (this->hasExactSourceBounds() && this->getSourceBounds(bounds))
	boundsExact = TRUE;
    else if (detached->getSourceBounds(bounds))
	boundsExact = detached->hasExactSourceBounds();
    else if (streamedComplete && this->getSourceBounds(bounds))
	boundsExact = this->hasExactSourceBounds();
    const std::shared_ptr<BObolCompactOccurrenceStream> nextStaging =
	stagedSourceStream && stagedSourceStream->stagedSourceByteCount() ?
	stagedSourceStream : std::shared_ptr<BObolCompactOccurrenceStream>();
    BObolCompactSourceProfile streamProfile;
    const bool haveStreamProfile = stagedSourceStream &&
	stagedSourceStream->getSourceProfile(streamProfile);
    const size_t streamExpectedCount = stagedSourceStream ?
	stagedSourceStream->getExpectedCount() : 0;
    if ((haveStreamProfile && !streamExpectedCount) ||
	(streamExpectedCount && !this->acceptsCompactStreamCertification(
	    streamExpectedCount, haveStreamProfile ? &streamProfile : nullptr)))
	return 0;
    std::vector<SoNode *> nextChildren;
    nextChildren.reserve(size_t(this->getNumChildren()));
    for (int i = 0; i < this->getNumChildren(); ++i) {
	SoNode *child = this->getChild(i);
	if (source_retains_realized_child(child))
	    nextChildren.push_back(child);
    }
    auto childRemoval = prepare_non_auxiliary_child_removal(this);
    PreparedSourceRealization publication(*this, bounds, boundsExact,
	!this->d->compactSelectedPaths.empty(), this->realizationRoleFlags.getValue());
    if (effects) {
	const std::vector<SoNode *> changedNodes;
	effects->prepare(*this, true, SbString(), nextChildren, changedNodes);
    }
    if (streamedComplete) {
	/* Once the complete source index has been adopted, all leaf coverage is
	 * available for the next retained assembly update.  Retire the temporary
	 * whole-target extent here as well as at presentation time: headless
	 * clients have no assembly callback, and leaving it visible there makes a
	 * completed draw permanently report a structural fallback. */
	std::vector<size_t> retiredOverviewEntries;
	retiredOverviewEntries.reserve(overviewCount);
	for (size_t i = 0; i < current->entries.size(); ++i) {
	    const BObolCompactInstanceEntry &entry = current->entries[i];
	    if ((entry.authoredVisible || entry.visible) &&
		BU_STR_EQUAL(entry.shapeSummary.recordRole.getString(),
		    "lod-overview"))
		retiredOverviewEntries.push_back(i);
	}
	current->hiddenInstances.reserve(current->hiddenInstances.size() +
	    retiredOverviewEntries.size());
	for (size_t i : retiredOverviewEntries) {
	    BObolCompactInstanceEntry &entry = current->entries[i];
	    if (entry.visible)
		current->hiddenInstances.push_back(entry.instance);
	    entry.authoredVisible = FALSE;
	    entry.visible = FALSE;
	    entry.visibilityRevision = compact_next_revision(entry.visibilityRevision);
	    compact_sync_shape_summary_state(entry);
	}
	if (!retiredOverviewEntries.empty()) {
	    this->markCompiledAssemblyDirty();
	    this->markCadBatchDirty(retiredOverviewEntries);
	}
	this->d->compactOverviewState =
	    BObolCompactOccurrenceRegistryState::OverviewState::Retired;
	/*
	 * The live compact registry is now the authoritative detached result.
	 * A source/input/database invalidation which started this realization
	 * may have revoked its request epoch while validated batches were being
	 * merged.  Successful adoption must publish a current epoch again;
	 * otherwise all request-bearing entries remain present but
	 * hasDisplayMeshLodRequests() reports false forever.  Wireframe then has
	 * no shaded fallback source and can never refine past the prefixes which
	 * happened to arrive before this handoff.
	 */
	const SbBool hadCurrentLodContract =
	    this->hasDisplayMeshLodRequests();
	this->d->displayMeshLodContractRevisionValid =
	    current->sourceMeshRequestCount > 0 ? TRUE : FALSE;
	this->d->displayMeshLodContractSourceRevision =
	    this->sourceRevision.getValue();
	this->d->displayMeshLodContractInputsRevision =
	    this->inputsRevision.getValue();
	if (!hadCurrentLodContract &&
	    this->d->displayMeshLodContractRevisionValid)
	    this->markDisplayMeshLodDirty();
	if (getenv("BOBOL_LOD_TRACE_SOURCE_CONTRACT"))
	    bu_log("BObol LoD source contract authoritative adoption path=%s "
		   "entries=%zu requests=%zu restored=%d\n",
		   this->path.getValue().getString(), current->entries.size(),
		   current->sourceMeshRequestCount,
		   (!hadCurrentLodContract &&
		    this->d->displayMeshLodContractRevisionValid) ? 1 : 0);

    } else {
	/*
	 * Camera revision is not source-content identity.  A user commonly moves
	 * the view while this worker replaces a structural frontier with retained
	 * native/PoP-backed occurrences.  Rejecting that valid result strands
	 * depth-capped combination boxes; the active PoP prefix is selected from
	 * the current view independently after adoption.
	 */
	std::unique_ptr<BObolCompactInstanceIndex> next(detached->d->compactIndex);
	detached->d->compactIndex = NULL;
	if (getenv("BOBOL_DRAW_TIMING")) {
	    const bool verbose = getenv("BOBOL_DRAW_TIMING_VERBOSE") != NULL;
	    size_t unresolved = 0;
	    for (const BObolCompactInstanceEntry &entry : next->entries) {
		if (!BU_STR_EQUAL(entry.shapeSummary.geometryKind.getString(),
			"aabb") &&
		    !BU_STR_EQUAL(entry.shapeSummary.geometryKind.getString(),
			"obb"))
		    continue;
		unresolved++;
		if (verbose)
		    bu_log("[obol-timing] deferred terminal proxy: path=%s "
			"mesh-request=%d lod-backed=%d\n",
			entry.semantic.path.getString(),
			entry.sourceMeshRequestValid ? 1 : 0,
			entry.lodBacked ? 1 : 0);
	    }
	    if (unresolved)
		bu_log("[obol-timing] deferred terminal proxies: n=%zu\n",
		    unresolved);
	}
	this->installCompactInstanceIndex(std::move(next), TRUE,
	    detached->d->compactHandleSourceId);
    }
    this->d->compactStagedSourceStream = nextStaging;
    if (streamExpectedCount)
	(void)this->certifyCompactStream(streamExpectedCount,
	    haveStreamProfile ? &streamProfile : nullptr);
    publication.commit();
    if (childRemoval)
	childRemoval->commit();
    if (effects)
	effects->commit(true);
    std::exception_ptr notificationFailure;
    publication.notify(notificationFailure);
    try {
	if (childRemoval)
	    childRemoval->notify();
    } catch (...) {
	if (!notificationFailure)
	    notificationFailure = std::current_exception();
    }
    if (effects)
	effects->notify(notificationFailure);
    if (notificationFailure)
	std::rethrow_exception(notificationFailure);
    if (getenv("BOBOL_DRAW_TIMING") && streamedComplete)
	bu_log("[obol-timing] deferred adoption preserved streamed index: n=%zu overview=%zu\n",
	    authoritativeCount, overviewCount);
    return authoritativeCount > static_cast<size_t>(INT_MAX) ? INT_MAX :
	static_cast<int>(authoritativeCount);
}

static SbString
retarget_record_path(const SbString &path, const char *oldSourcePath, const char *newSourcePath)
{
    if (path.getLength() && BU_STR_EQUAL(oldSourcePath, newSourcePath)) return path;
    const std::string_view record(database_source_skip_leading_slash(path.getString()));
    const std::string_view previous(database_source_skip_leading_slash(oldSourcePath));
    if (record.empty() || record == previous) return newSourcePath;
    if (!previous.empty() && record.size() > previous.size() &&
	record.substr(0, previous.size()) == previous && record[previous.size()] == '/') {
	SbString next(newSourcePath);
	next += record.data() + previous.size();
	return next;
    }
    return path;
}

template <typename ShapeT>
static void
retarget_realized_shape_source(ShapeT *shape,
			       const SoBRLDatabaseSource *source,
			       const char *oldSourcePath,
			       const char *newSourcePath,
			       uint32_t revision)
{
    if (!shape || !newSourcePath)
	return;

    const std::string oldNameStorage =
	stable_name_from_path(oldSourcePath, 0);
    const char *oldName = oldNameStorage.empty() ? NULL :
			  oldNameStorage.c_str();
    const std::string newNameStorage =
	stable_name_from_path(newSourcePath, 1);
    const char *newName = newNameStorage.empty() ? newSourcePath :
			  newNameStorage.c_str();

    const std::string recordRole = shape->recordRole.getValue().getString();
    const int auxiliary = BU_STR_EQUAL(recordRole.c_str(), "auxiliary");
    const std::string shapeSourceName =
	shape->sourceName.getValue().getString();
    const std::string displayName = shape->displayName.getValue().getString();
    const std::string geometryName =
	shape->geometryName.getValue().getString();
    const std::string sourceDisplayName = source ?
					  source->displayName.getValue().getString() : "";

    const SbString nextPath = retarget_record_path(shape->sourcePath.getValue(), oldSourcePath, newSourcePath);
    database_source_assign_string(shape->sourcePath, nextPath);
    shape->sourceId = revision;
    if (!auxiliary &&
	(shapeSourceName.empty() ||
	 (oldName && BU_STR_EQUAL(shapeSourceName.c_str(), oldName)) ||
	 (oldSourcePath && BU_STR_EQUAL(shapeSourceName.c_str(), oldSourcePath))))
	database_source_assign_string(shape->sourceName, newName);
    if (!auxiliary && !sourceDisplayName.empty())
	database_source_assign_string(shape->displayName,
				      sourceDisplayName.c_str());
    else if (!auxiliary &&
	     (displayName.empty() ||
	      (oldName && BU_STR_EQUAL(displayName.c_str(), oldName)) ||
	      (oldSourcePath && BU_STR_EQUAL(displayName.c_str(), oldSourcePath))))
	database_source_assign_string(shape->displayName, newName);
    if (!auxiliary &&
	(geometryName.empty() ||
	 (oldName && BU_STR_EQUAL(geometryName.c_str(), oldName)) ||
	 (oldSourcePath && BU_STR_EQUAL(geometryName.c_str(), oldSourcePath))))
	database_source_assign_string(shape->geometryName, newName);

    SbString identity = source_record_identity(source, nextPath.getString());
    if (auxiliary) {
	const char *auxName = shape->geometryName.getValue().getString();
	if (auxName && auxName[0]) {
	    identity += "::";
	    identity += auxName;
	}
    }
    database_source_assign_string(shape->sourceIdentity, identity);
    database_source_assign_string(shape->cacheIdentity,
				  record_identity_with_revision(identity.getString(), revision));
    sync_shape_owner_state(shape, source);
}

static void
retarget_material_object_source(SoBRLMaterialObject *object,
				const char *oldSourcePath,
				const char *newSourcePath,
				uint32_t revision)
{
    if (!object || !newSourcePath)
	return;

    const std::string oldNameStorage =
	stable_name_from_path(oldSourcePath, 0);
    const char *oldName = oldNameStorage.empty() ? NULL :
			  oldNameStorage.c_str();
    const std::string newNameStorage =
	stable_name_from_path(newSourcePath, 1);
    const char *newName = newNameStorage.empty() ? newSourcePath :
			  newNameStorage.c_str();

    const char *sourceName = object->sourceName.getValue().getString();
    object->sourcePath = retarget_record_path(object->sourcePath.getValue(), oldSourcePath, newSourcePath);
    object->sourceId = revision;
    if (!sourceName || !sourceName[0] ||
	(oldName && BU_STR_EQUAL(sourceName, oldName)) ||
	(oldSourcePath && BU_STR_EQUAL(sourceName, oldSourcePath)))
	object->sourceName = newName;
}

int
SoBRLDatabaseSource::retargetDatabaseSource(const char *sourcePath,
	uint32_t revision)
{
    return this->retargetDatabaseSourceInstance(sourcePath, sourcePath,
	    revision);
}

SbBool
SoBRLDatabaseSource::needsRealization(void) const
{
    return this->stale.getValue() ||
	   this->realizedSourceRevision.getValue() != this->sourceRevision.getValue() ||
	   this->realizedInputsRevision.getValue() != this->inputsRevision.getValue() ||
	   this->realizedViewRevision.getValue() != this->viewRevision.getValue();
}

static int
source_has_view_lod_payload(SoAction *action, SoNode *node)
{
    if (!action || !node)
	return 0;

    /* A compact source is rendered and picked from its occurrence registry.
     * Per-shape payloads may still be present in the view state while the LoD
     * controller transitions to aggregate payloads, but traversing those
     * shapes would discard the compact representation and recreate the large
     * Coin render graph it is intended to replace.  A source-level payload is
     * handled by cad_view_lod_assembly_for_action before this query. */
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoBRLDatabaseSource *source =
	    static_cast<SoBRLDatabaseSource *>(node);
	if (source->hasCompactInstanceIndex())
	    return 0;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *mesh = static_cast<SoBRLMeshShape *>(node);
	return bobol_view_lod_mesh_for_action(action, mesh) ||
	       bobol_view_lod_proxy_for_action(action, mesh);
    }
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return 0;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++) {
	if (source_has_view_lod_payload(action, group->getChild(i)))
	    return 1;
    }
    return 0;
}

void
SoBRLDatabaseSource::GLRender(SoGLRenderAction *action)
{
    if (bobol_cad_batch_source_suppressed(this))
	return;
    if (SoBRLCadAssembly *viewCad =
	    cad_view_lod_assembly_for_action(action, this)) {
	viewCad->render(action);
	return;
    }

    if (source_has_view_lod_payload(action, this)) {
	inherited::GLRender(action);
	return;
    }

    if (this->syncCompiledAssembly() && this->d->compiledAssembly) {
	this->d->compiledAssembly->render(action);
	return;
    }

    inherited::GLRender(action);
}

void
SoBRLDatabaseSource::GLRenderBelowPath(SoGLRenderAction *action)
{
    if (bobol_cad_batch_source_suppressed(this))
	return;
    if (SoBRLCadAssembly *viewCad =
	    cad_view_lod_assembly_for_action(action, this)) {
	viewCad->render(action);
	return;
    }

    if (source_has_view_lod_payload(action, this)) {
	inherited::GLRenderBelowPath(action);
	return;
    }

    if (this->syncCompiledAssembly() && this->d->compiledAssembly) {
	this->d->compiledAssembly->render(action);
	return;
    }

    inherited::GLRenderBelowPath(action);
}

void
SoBRLDatabaseSource::callback(SoCallbackAction *action)
{
    /* Non-GL renderers traverse the retained aggregate through its one node. */
    if (this->hasCompactInstanceIndex())
	(void)this->syncCompiledAssembly();
    inherited::callback(action);
}

void
SoBRLDatabaseSource::getBoundingBox(SoGetBoundingBoxAction *action)
{
    if (SoBRLCadAssembly *viewCad =
	    cad_view_lod_assembly_for_action(action, this)) {
	viewCad->getBounds(action);
	return;
    }

    if (this->hasCompactInstanceIndex()) {
	SbBox3f bounds;
	bounds.makeEmpty();
	if (this->visible.getValue()) {
	    for (const BObolCompactInstanceEntry &entry :
		 this->d->compactIndex->entries) {
		if (!entry.visible || !entry.geometry)
		    continue;
		bounds.extendBy(database_source_transform_bounds(
		    compact_part_geometry_bounds(entry.geometry),
		    entry.localToSource));
	    }
	}
	if (!bounds.isEmpty()) {
	    action->extendBy(bounds);
	    action->setCenter(bounds.getCenter(), TRUE);
	}
	return;
    }

    inherited::getBoundingBox(action);
}

void
SoBRLDatabaseSource::rayPick(SoRayPickAction *action)
{
    if (SoBRLCadAssembly *viewCad =
	    cad_view_lod_assembly_for_action(action, this)) {
	viewCad->pickRay(action);
	return;
    }

    if (this->hasCompactInstanceIndex() &&
	this->syncCompiledAssembly() && this->d->compiledAssembly) {
	this->d->compiledAssembly->pickRay(action);
	return;
    }

    inherited::rayPick(action);
}

SbBool
SoBRLDatabaseSource::realizeDatabaseWireframe(void)
{
    return this->realizeDatabaseWireframe(NULL);
}
SbBool
SoBRLDatabaseSource::realizeDatabaseWireframe(
    BObolCompactOccurrenceStream *stream)
{
    BObolDatabaseSourceRealizationCache cache;
    if (source_uses_evaluated_wire_realization(this))
	return bobol_database_source_realize_wireframe_with_cache(this, &cache);
    return bobol_database_source_realize_wireframe_compact_with_cache(
	this, &cache, stream) > 0 ? TRUE : FALSE;
}

static void
mark_source_realized_current(SoBRLDatabaseSource *source)
{
    if (!source)
	return;

    source->realizedRevision = source->sourceRevision.getValue();
    source->realizedSourceRevision = source->sourceRevision.getValue();
    source->realizedInputsRevision = source->inputsRevision.getValue();
    source->realizedViewRevision = source->viewRevision.getValue();
    source->realizationStatus = SoBRLDatabaseSource::REALIZED;
    source->realizationDiagnostic = "";
    source->realizationIdentity = source_realization_identity(source);
    source->stale = FALSE;
    source->staleReason = SoBRLDatabaseSource::STALE_NONE;
}

int
bobol_database_source_construct_wireframe_compact_with_cache(
    SoBRLDatabaseSource *source,
    BObolDatabaseSourceRealizationCache *cache,
    BObolCompactOccurrenceStream *stream)
{
    BObolPerformanceTimer timer(BOBOL_PERF_WIRE_REALIZE_US);
    if (timer.active())
	bobol_performance_counter_add(BOBOL_PERF_WIRE_REALIZE_CALLS, 1);

    BObolDatabaseSourceRealizationCache localCache;
    if (!cache)
	cache = &localCache;

    if (!source)
	return -1;

    if (source_uses_evaluated_path_realization(source))
	return 0;

    source->d->compactHandleSourceId = source_stable_compact_handle_id(source);

    source->realizationDiagnostic = "";
    if (!source->d->dbip) {
	source->realizationDiagnostic = "database source has no database";
	return -1;
    }

    std::string treeNameStorage =
	database_lookup_path_from_source_path(source->path.getValue());
    const char *treeName = treeNameStorage.c_str();
    if (!treeName[0]) {
	source->realizationDiagnostic = "database source path is empty";
	return -1;
    }

    (void)treeName;
    const uint32_t revision = source->sourceRevision.getValue();
    int directRealized = 0;
    {
	BObolPerformanceTimer directTimer(BOBOL_PERF_DIRECT_LEAF_US);
	if (directTimer.active())
	    bobol_performance_counter_add(BOBOL_PERF_DIRECT_LEAF_CALLS, 1);
	directRealized = realize_direct_leaf_wireframe_compact(source, cache,
	    revision, stream);
	if (directTimer.active()) {
	    if (directRealized > 0) {
		bobol_performance_counter_add(
		    BOBOL_PERF_DIRECT_LEAF_REALIZED, 1);
	    } else if (directRealized < 0) {
		bobol_performance_counter_add(
		    BOBOL_PERF_DIRECT_LEAF_FAILED, 1);
	    } else {
		bobol_performance_counter_add(
		    BOBOL_PERF_DIRECT_LEAF_FALLBACK, 1);
	    }
	}
    }

    if (directRealized > 0) {
	SbBox3f semanticBounds;
	if (source_bounds_from_database_path(source, treeName, semanticBounds))
	    (void)source->setSourceBoundsState(TRUE, semanticBounds.getMin(),
		semanticBounds.getMax(), TRUE);
	mark_source_realized_current(source);
	return 1;
    }
    if (directRealized < 0) {
	remove_non_auxiliary_children(source);
	source->discardCompactInstanceHistory();
	source->realizationIdentity = "";
	return -1;
    }

    struct db_tree_state init_state;
    db_init_db_tree_state(&init_state, source->d->dbip);
    init_state.ts_stop_at_regions = 0;

    BObolMaterialColorSweep materialSweep(source->d->dbip);
    realize_walk_data data;
    data.source = source;
    data.cache = cache;
    data.revision = revision;
    data.compact_index = std::make_unique<BObolCompactInstanceIndex>();
    data.stream_sink = stream;
    data.material_sweep = &materialSweep;

    const char *av[1] = {treeName};
    const int ret = db_walk_tree_leaf_instances(source->d->dbip, 1, av, 1,
	&init_state, NULL, NULL, realize_leaf, &data);
    db_free_db_tree_state(&init_state);

    if (ret < 0 || data.realized_shapes <= 0 || data.failed_shapes > 0 ||
	data.compact_unsupported || data.compact_index->entries.empty()) {
	const size_t compactEntryCount = data.compact_index->entries.size();
	remove_non_auxiliary_children(source);
	source->discardCompactInstanceHistory();
	source->realizationIdentity = "";
	if (data.diagnostic.getLength() > 0)
	    source->realizationDiagnostic = data.diagnostic;
	else {
	    SbString diagnostic;
	    diagnostic.sprintf(
		"%s: compact wireframe realization produced no usable "
		"occurrences (walk=%d, leaves=%zu, realized=%d, failed=%d, "
		"unsupported=%d, entries=%zu, instance_ids=%d)",
		treeName, ret, data.visited_leaves, data.realized_shapes,
		data.failed_shapes, data.compact_unsupported,
		compactEntryCount,
		db_comb_instance_ids_get(source->d->dbip));
	    source->realizationDiagnostic = diagnostic;
	}
	return ret < 0 || data.failed_shapes > 0 ? -1 : 0;
    }

    source->installCompactInstanceIndex(std::move(data.compact_index), TRUE,
	source->d->compactHandleSourceId);
    source->markCompiledAssemblyDirty();

    if (stream)
	(void)close_compact_stream_contract(source, stream);

    const SbBool preserveExact = source->hasExactSourceBounds();

    SbBox3f semanticBounds;
    if (!preserveExact && source_bounds_from_database_path(source, treeName,
	    semanticBounds))
	(void)source->setSourceBoundsState(TRUE, semanticBounds.getMin(),
	    semanticBounds.getMax(), TRUE);
    else if (data.compact_bounds_valid && !data.compact_bounds.isEmpty() &&
	!preserveExact)
	set_source_bounds_from_presentation(source, data.compact_bounds);
    else if (!preserveExact &&
	(!data.compact_bounds_valid || data.compact_bounds.isEmpty()))
	source->clearSourceBounds();
    mark_source_realized_current(source);
    bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_SOURCES, 1);
    bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_INSTANCES,
	static_cast<uint64_t>(source->d->compactIndex->entries.size()));
    return 1;
}

SbBool
bobol_database_source_construct_wireframe_with_cache(
    SoBRLDatabaseSource *source,
    BObolDatabaseSourceRealizationCache *cache)
{
    BObolPerformanceTimer timer(BOBOL_PERF_WIRE_REALIZE_US);
    if (timer.active())
	bobol_performance_counter_add(BOBOL_PERF_WIRE_REALIZE_CALLS, 1);

    BObolDatabaseSourceRealizationCache localCache;
    if (!cache)
	cache = &localCache;

    if (!source)
	return FALSE;

    source->realizationDiagnostic = "";
    if (!source->d->dbip) {
	source->realizationDiagnostic = "database source has no database";
	return FALSE;
    }

    std::string treeNameStorage =
	database_lookup_path_from_source_path(source->path.getValue());
    const char *treeName = treeNameStorage.c_str();
    if (!treeName[0]) {
	source->realizationDiagnostic = "database source path is empty";
	return FALSE;
    }

    remove_non_auxiliary_children(source);
    source->discardCompactInstanceHistory();
    (void)remove_source_placement_transform(source);

    const uint32_t revision = source->sourceRevision.getValue();
    if (source_uses_evaluated_wire_realization(source)) {
	if (realize_evaluated_wire_source(source, revision) > 0) {
	    mark_source_realized_current(source);
	    /* Evaluated wire is a terminal representation, not a construction
	     * geometry union.  Prefer librt's representation-independent path
	     * bound so progressive and explicit autoview share one camera.  Some
	     * valid evaluated paths (for example air-only assemblies) are omitted
	     * by rt_obj_bounds' ordinary solid policy; in that case the evaluated
	     * result itself is the authoritative representation and its bound can
	     * safely close the readiness contract. */
	    SbBox3f semanticBounds;
	    if (source_bounds_from_database_path(source, treeName,
		    semanticBounds))
		(void)source->setSourceBoundsState(TRUE,
		    semanticBounds.getMin(), semanticBounds.getMax(), TRUE);
	    else
		update_source_bounds_from_realized_geometry(source, TRUE);
	    source->syncRealizedShapeOwnerState();
	    return TRUE;
	}
	remove_non_auxiliary_children(source);
	source->realizationIdentity = "";
	return FALSE;
    }

    int directRealized = 0;
    {
	BObolPerformanceTimer directTimer(BOBOL_PERF_DIRECT_LEAF_US);
	if (directTimer.active())
	    bobol_performance_counter_add(BOBOL_PERF_DIRECT_LEAF_CALLS, 1);
	directRealized = realize_direct_leaf_wireframe(source, cache, revision);
	if (directTimer.active()) {
	    if (directRealized > 0) {
		bobol_performance_counter_add(
		    BOBOL_PERF_DIRECT_LEAF_REALIZED, 1);
	    } else if (directRealized < 0) {
		bobol_performance_counter_add(
		    BOBOL_PERF_DIRECT_LEAF_FAILED, 1);
	    } else {
		bobol_performance_counter_add(
		    BOBOL_PERF_DIRECT_LEAF_FALLBACK, 1);
	    }
	}
    }
    if (directRealized > 0) {
	SbBox3f semanticBounds;
	if (source_bounds_from_database_path(source, treeName, semanticBounds))
	    (void)source->setSourceBoundsState(TRUE, semanticBounds.getMin(),
		semanticBounds.getMax(), TRUE);
	mark_source_realized_current(source);
	SbBox3f sourceBounds;
	if (!source->getSourceBounds(sourceBounds))
	    update_source_bounds_from_realized_geometry(source);
	source->syncRealizedShapeOwnerState();
	return TRUE;
    }
    if (directRealized < 0) {
	remove_non_auxiliary_children(source);
	source->realizationIdentity = "";
	return FALSE;
    }

    struct db_tree_state init_state;
    db_init_db_tree_state(&init_state, source->d->dbip);
    init_state.ts_stop_at_regions = 0;

    struct realize_walk_data data;
    data.source = source;
    data.cache = cache;
    data.revision = source->sourceRevision.getValue();
    data.realized_shapes = 0;
    data.failed_shapes = 0;

    const char *av[1] = { treeName };
    int ret = db_walk_tree_leaf_instances(source->d->dbip, 1, av, 1, &init_state,
					  NULL, NULL, realize_leaf, &data);
    db_free_db_tree_state(&init_state);

    if (ret < 0 || data.realized_shapes <= 0 || data.failed_shapes > 0) {
	remove_non_auxiliary_children(source);
	source->realizationIdentity = "";
	if (data.diagnostic.getLength() > 0) {
	    source->realizationDiagnostic = data.diagnostic;
	} else if (data.realized_shapes <= 0) {
	    SbString msg;
	    msg.sprintf("%s: no drawable wireframe geometry realized", treeName);
	    source->realizationDiagnostic = msg;
	} else {
	    SbString msg;
	    msg.sprintf("%s: wireframe realization failed", treeName);
	    source->realizationDiagnostic = msg;
	}
	return FALSE;
    }

    mark_source_realized_current(source);
    update_source_bounds_from_realized_geometry(source);
    source->syncRealizedShapeOwnerState();
    return TRUE;
}

SbBool
SoBRLDatabaseSource::realizeDatabaseMesh(void)
{
    return this->realizeDatabaseMesh(NULL);
}

SbBool
SoBRLDatabaseSource::realizeDatabaseMesh(BObolCompactOccurrenceStream *stream)
{
    BObolDatabaseSourceRealizationCache cache;
    if (source_uses_evaluated_points_realization(this))
	return bobol_database_source_realize_mesh_with_cache(this, &cache);
    return bobol_database_source_realize_mesh_compact_with_cache(
	this, &cache, stream) > 0 ? TRUE : FALSE;
}

namespace {

struct compact_coverage_occurrence {
    SbMatrix localTransform = SbMatrix::identity();
    BObolRealizedShapeSummary summary;
    std::string assetPath;
    uint32_t occurrenceIndex = 0;
    int booleanOperation = SoBRLDatabaseSource::BOOLEAN_UNION;
};

struct compact_coverage_asset {
    struct directory *dp = NULL;
    std::string cacheKey;
    std::string assetPath;
    size_t estimatedWorkingSetBytes = 1;
    SbBox3f coverageBounds;
    std::shared_ptr<const Obol::PartGeometry> coverageGeometry;
    std::shared_ptr<const Obol::PartGeometry> geometry;
    SbMatrix proxyGeometryTransform = SbMatrix::identity();
    BObolSourceMeshRequest sourceMeshRequest;
    uint64_t sampleFingerprint = 0;
    bool sampleFingerprintValid = false;
    size_t vertexCount = 0;
    size_t faceCount = 0;
    unsigned char mode = 0;
    unsigned char orientation = 0;
    unsigned char flags = 0;
    bool coverageReady = false;
    bool lodEligible = false;
    /* True only when geometry is the standing proxy for a source-backed
     * triangle PoP asset.  BREP wire drawing has its own immutable
     * progressive line representation and must not also submit the shaded
     * tessellation as a wire overlay. */
    bool sourceMeshReady = false;
    /* A persisted canonical mapping is already an exact-reuse proof for this
     * asset.  It may supply a warm source contract immediately, but it must
     * not become the representative of a new cold proof group: its request
     * transform maps the canonical asset to this object, not this object to
     * another candidate.  Treating it as a representative without composing
     * those transforms can write a non-identity self mapping for the
     * canonical object and draw later warm meshes at the wrong location. */
    bool sourceMeshMappingCached = false;
    /* A syntactically valid BoT with no faces has no drawable occurrence.
     * Distinguish that terminal no-op from malformed nonempty data: the
     * former is excluded from the producer's final drawable census, while
     * the latter must still prevent authoritative stream adoption. */
    bool terminalEmpty = false;
    std::string realizedSourceType;
    std::string realizedGeometryKind;
    bool viewDependentCsgGeometry = false;
    /* A later asset with the same cheap rigid-invariant signature must not
     * launch its own potentially enormous PoP build until exact serialized
     * comparison has either proved reuse or rejected it. */
    bool deferSourceMeshContract = false;
    std::vector<compact_coverage_occurrence> deferredOccurrences;
    bool ready = false;
    std::once_flag coverageOnce;
    std::once_flag realizeOnce;
};

struct compact_coverage_work_item {
    compact_coverage_asset *asset = NULL;
    compact_coverage_occurrence occurrence;
};

/* Amortize the producer/consumer mutex without delaying useful feedback by
 * more than a tiny fraction of one frame.  Discovery order is not a rendering
 * contract; semantic identity is assigned before an item enters this queue. */
static const size_t compact_coverage_work_batch_size = 16;

struct compact_coverage_collect {
    SoBRLDatabaseSource *source = NULL;
    BObolDatabaseSourceRealizationCache *cache = NULL;
    BObolCompactOccurrenceStream *stream = NULL;
    BObolMaterialColorSweep *materialSweep = NULL;
    uint32_t revision = 0;
    std::vector<std::unique_ptr<compact_coverage_asset>> assets;
    /* A directory entry is the canonical object identity for the lifetime of
     * this detached realization database.  Hashing the complete persistent
     * cache key again for every occurrence duplicated librt's name lookup and
     * string work on the serial discovery path.  Assets retain cacheKey for
     * persistence; only this walk-local index uses the stable pointer. */
    std::unordered_map<struct directory *, size_t> assetIndices;
    std::unordered_set<std::string> seenInstances;
    std::unordered_map<std::string, uint32_t> occurrenceCounts;
    size_t occurrenceCount = 0;
    /* A complete union of leaf boxes is an exact source bound only when each
     * leaf bbox is itself tight and Boolean evaluation cannot shrink it.
     * BRep GetBBox may include untrimmed surface extent, while subtract and
     * intersect operations may remove an extremum.  Keep publishing the
     * conservative aggregate for immediate visual feedback, but require a
     * semantic librt bound before releasing the autoview/cache contract. */
    bool aggregateBoundsExact = true;
    std::mutex workMutex;
    std::condition_variable workReady;
    std::deque<compact_coverage_work_item> work;
    /* Written only by the hierarchy-walk producer, then transferred as one
     * bounded queue operation. */
    std::vector<compact_coverage_work_item> producerWork;
    /* Coverage is the latency-critical phase.  Retain its completed work
     * records here, then import full BoTs only after every leaf box and the
     * exact target extent have been published. */
    std::deque<compact_coverage_work_item> detailWork;
    /* BOT copies and BREP providers have their own thread-safe paths.  The
     * remaining primitive conversion stack may enter NMG/plot code with
     * process-global error state, so only that small portion is serialized
     * while imports and publication remain parallel. */
    std::mutex primitiveGeometryMutex;
    std::mutex reuseMutex;
    std::unordered_map<std::string,
	std::vector<compact_coverage_asset *>> reuseGroups;
    bool producerDone = false;
};

static size_t
compact_coverage_working_set_estimate(const struct db_i *dbip,
	const struct directory *dp)
{
    const bool brep = dp &&
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP;
    const size_t fixedBytes = brep ?
	64ULL * 1024ULL * 1024ULL : 8ULL * 1024ULL * 1024ULL;
    if (!dp)
	return fixedBytes;
    size_t encodedBytes = dp->d_len;
    if (dbip && db_version(dbip) < 5) {
	if (encodedBytes > SIZE_MAX / sizeof(union record))
	    return SIZE_MAX;
	encodedBytes *= sizeof(union record);
    }
    /*
     * A BoT coverage worker transiently owns the database read buffer, the
     * decoded rt_bot_internal arrays, and modest import bookkeeping.  The
     * previous 16x multiplier charged Lucy as 10.8 GiB for roughly 675 MiB of
     * source arrays.  Since the global governor permits an oversized task
     * only when it runs alone, that estimate also serialized ordinary
     * tens-of-megabytes vehicle parts despite ample RAM and CPUs.
     *
     * Three serialized copies plus fixed overhead conservatively covers the
     * encoded record, decoded arrays, and decode scratch for current v5 BoTs.
     * The retained source arrays then move into the separately bounded staged
     * lease window; downstream PoP construction takes its own face/point
     * reservation before it can overlap another import.
     */
    const size_t copyFactor = brep ? 8 : 3;
    if (encodedBytes > (SIZE_MAX - fixedBytes) / copyFactor)
	return SIZE_MAX;
    return encodedBytes * copyFactor + fixedBytes;
}

struct compact_coverage_vertex_bounds {
    point_t minimum = {INFINITY, INFINITY, INFINITY};
    point_t maximum = {-INFINITY, -INFINITY, -INFINITY};
    size_t finiteVertices = 0;
};

/* Keep nested large-mesh scans from multiplying the ordinary coverage worker
 * pool.  A small collection of large BoTs still uses otherwise-idle CPUs,
 * while a 50k/150k collection of ordinary leaves retains outer parallelism. */
static constexpr size_t compact_coverage_parallel_vertex_threshold =
    1024 * 1024;
static constexpr size_t compact_coverage_parallel_vertex_workers = 4;
static constexpr size_t compact_coverage_parallel_vertex_scans = 2;
static std::atomic<size_t> compact_coverage_active_vertex_scans(0);

static compact_coverage_vertex_bounds
compact_coverage_scan_serialized_vertices(const unsigned char *verticesBody,
	 size_t firstVertex, size_t vertexCount, size_t vertexStride)
{
    compact_coverage_vertex_bounds result;
    if (!verticesBody || !vertexCount)
	return result;

    static constexpr size_t vertexBlock = 256;
    double decoded[vertexBlock * ELEMENTS_PER_POINT];
    const unsigned char *point = verticesBody + firstVertex * vertexStride;
    for (size_t first = 0; first < vertexCount; first += vertexBlock) {
	const size_t count = std::min(vertexBlock, vertexCount - first);
	bu_cv_ntohd(reinterpret_cast<unsigned char *>(decoded), point,
	    count * ELEMENTS_PER_POINT);
	point += count * vertexStride;
	for (size_t i = 0; i < count; ++i) {
	    const double *value = decoded + i * ELEMENTS_PER_POINT;
	    if (!std::isfinite(value[X]) || !std::isfinite(value[Y]) ||
		!std::isfinite(value[Z]))
		continue;
	    for (int axis = 0; axis < 3; ++axis) {
		result.minimum[axis] = std::min(result.minimum[axis],
		    value[axis]);
		result.maximum[axis] = std::max(result.maximum[axis],
		    value[axis]);
	    }
	    result.finiteVertices++;
	}
    }
    return result;
}

static void
compact_coverage_merge_vertex_bounds(compact_coverage_vertex_bounds &result,
	const compact_coverage_vertex_bounds &candidate)
{
    if (!candidate.finiteVertices)
	return;
    if (!result.finiteVertices) {
	result = candidate;
	return;
    }
    for (int axis = 0; axis < 3; ++axis) {
	result.minimum[axis] = std::min(result.minimum[axis],
	    candidate.minimum[axis]);
	result.maximum[axis] = std::max(result.maximum[axis],
	    candidate.maximum[axis]);
    }
    result.finiteVertices += candidate.finiteVertices;
}

static bool
compact_coverage_try_acquire_parallel_vertex_scan(void)
{
    size_t active = compact_coverage_active_vertex_scans.load();
    while (active < compact_coverage_parallel_vertex_scans) {
	if (compact_coverage_active_vertex_scans.compare_exchange_weak(active,
		active + 1))
	    return true;
    }
    return false;
}

/*
 * Decode only the fixed BoT header and vertex array needed for an AABB.
 *
 * rt_db_get_internal necessarily allocates and converts faces, thickness,
 * normals, UVs, and their index arrays as well.  That is appropriate for the
 * later PoP source contract, but it made first visual coverage of a many-part
 * database wait behind essentially the whole mesh import.  V5 keeps vertices
 * first in the BoT body, so this bounded decoder can publish every leaf box
 * and the exact draw-target extent before the detail phase pays those costs.
 *
 * A detached read-only realization database normally supplies a borrowed
 * memory-mapped record, so this pass scans the serialized vertices in place
 * without allocating or copying the faces which follow them.  Non-mapped v5
 * databases retain the db_get_external fallback.  V4 retains the complete
 * import fallback because its record layout does not offer the same compact
 * body contract.
 */
enum class CompactCoverageBotStatus : uint8_t {
    INVALID = 0,
    EMPTY,
    READY
};

static CompactCoverageBotStatus
compact_coverage_bot_bounds(struct db_i *dbip, struct directory *dp,
	SbBox3f &bounds, size_t &vertexCount, size_t &faceCount,
	unsigned char &mode, unsigned char &orientation, unsigned char &flags,
	uint64_t &sampleFingerprint, bool &sampleFingerprintValid,
	const struct bu_mapped_file *mappedFile = NULL)
{
    bounds.makeEmpty();
    vertexCount = 0;
    faceCount = 0;
    mode = 0;
    orientation = 0;
    flags = 0;
    sampleFingerprint = 0;
    sampleFingerprintValid = false;
    if (!dbip || !dp ||
	dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT)
	return CompactCoverageBotStatus::INVALID;

    if (db_version(dbip) != 5) {
	struct rt_db_internal intern;
	RT_DB_INTERNAL_INIT(&intern);
	if (rt_db_get_internal(&intern, dp, dbip, NULL) < 0)
	    return CompactCoverageBotStatus::INVALID;
	const struct rt_bot_internal *bot =
	    intern.idb_type == ID_BOT && intern.idb_ptr ?
	    static_cast<const struct rt_bot_internal *>(intern.idb_ptr) :
	    NULL;
	CompactCoverageBotStatus status = CompactCoverageBotStatus::INVALID;
	if (bot) {
	    vertexCount = bot->num_vertices;
	    faceCount = bot->num_faces;
	    mode = bot->mode;
	    orientation = bot->orientation;
	    flags = bot->bot_flags;
	    if (!faceCount) {
		status = CompactCoverageBotStatus::EMPTY;
	    } else if (vertexCount) {
		BObolSourceMeshRequest request;
		if (cad_source_mesh_request_from_bot(request, bot)) {
		    bounds = request.bounds;
		    sampleFingerprintValid =
			compact_stream_lod_sample_fingerprint(
			    sampleFingerprint, bot);
		    status = bounds.isEmpty() ?
			CompactCoverageBotStatus::INVALID :
			CompactCoverageBotStatus::READY;
		}
	    }
	}
	rt_db_free_internal(&intern);
	return status;
    }

    struct bu_external external;
    BU_EXTERNAL_INIT(&external);
    size_t serializedBytes = 0;
    const unsigned char *serialized = db_external_view(
	dbip, dp, &serializedBytes);
    bool borrowed = serialized != NULL;
    if (!serialized && mappedFile && mappedFile->buf && dp->d_addr >= 0 &&
	static_cast<size_t>(dp->d_addr) <= mappedFile->buflen &&
	dp->d_len <= mappedFile->buflen - static_cast<size_t>(dp->d_addr)) {
	serialized = static_cast<const unsigned char *>(mappedFile->buf) +
	    static_cast<size_t>(dp->d_addr);
	serializedBytes = dp->d_len;
	borrowed = true;
    }
    if (!borrowed) {
	if (db_get_external(&external, dp, dbip) < 0)
	    return CompactCoverageBotStatus::INVALID;
	serialized = external.ext_buf;
	serializedBytes = external.ext_nbytes;
    }

    const size_t headerBytes = sizeof(struct db5_ondisk_header);
    size_t encodedObjectUnits = 0;
    size_t objectLengthBytes = 0;
    bool serializedEnvelopeValid = serializedBytes >= headerBytes;
    if (serializedEnvelopeValid) {
	const int width = (serialized[1] &
	    DB5HDR_HFLAGS_OBJECT_WIDTH_MASK) >>
	    DB5HDR_HFLAGS_OBJECT_WIDTH_SHIFT;
	const size_t widthBytes = static_cast<size_t>(1) << width;
	serializedEnvelopeValid = width >= DB5HDR_WIDTHCODE_8BIT &&
	    width <= DB5HDR_WIDTHCODE_64BIT &&
	    headerBytes <= serializedBytes &&
	    widthBytes <= serializedBytes - headerBytes;
	if (serializedEnvelopeValid) {
	    (void)db5_decode_length(&encodedObjectUnits,
		serialized + headerBytes, width);
	    serializedEnvelopeValid =
		encodedObjectUnits <= SIZE_MAX / 8;
	    objectLengthBytes = serializedEnvelopeValid ?
		encodedObjectUnits * 8 : 0;
	    serializedEnvelopeValid = objectLengthBytes >= headerBytes &&
		objectLengthBytes <= serializedBytes;
	}
    }

    struct db5_raw_internal raw;
    const bool rawValid =
	serializedEnvelopeValid &&
	db5_get_raw_internal_ptr(&raw, serialized) != NULL &&
	raw.object_length == objectLengthBytes &&
	raw.major_type == DB5_MAJORTYPE_BRLCAD &&
	raw.minor_type == DB5_MINORTYPE_BRLCAD_BOT &&
	raw.body.ext_buf && raw.body.ext_nbytes >=
	    2 * SIZEOF_NETWORK_LONG + 3 &&
	raw.body.ext_buf >= serialized &&
	static_cast<size_t>(raw.body.ext_buf - serialized) <
	    raw.object_length &&
	raw.body.ext_nbytes <= raw.object_length -
	    static_cast<size_t>(raw.body.ext_buf - serialized) - 1;
    if (!rawValid) {
	if (!borrowed)
	    bu_free_external(&external);
	return CompactCoverageBotStatus::INVALID;
    }

    const unsigned char *body = raw.body.ext_buf;
    const size_t vertices = static_cast<size_t>(BU_GLONG(body));
    const size_t faces = static_cast<size_t>(
	BU_GLONG(body + SIZEOF_NETWORK_LONG));
    vertexCount = vertices;
    faceCount = faces;
    const size_t fixedBytes = 2 * SIZEOF_NETWORK_LONG + 3;
    const size_t vertexStride =
	SIZEOF_NETWORK_DOUBLE * ELEMENTS_PER_POINT;
    const size_t faceStride = 3 * SIZEOF_NETWORK_LONG;
    if (!faces) {
	if (!borrowed)
	    bu_free_external(&external);
	return CompactCoverageBotStatus::EMPTY;
    }
    if (!vertices || vertices > (SIZE_MAX - fixedBytes) / vertexStride ||
	fixedBytes + vertices * vertexStride > raw.body.ext_nbytes ||
	faces > (raw.body.ext_nbytes -
	    (fixedBytes + vertices * vertexStride)) / faceStride) {
	if (!borrowed)
	    bu_free_external(&external);
	return CompactCoverageBotStatus::INVALID;
    }

    orientation = body[2 * SIZEOF_NETWORK_LONG];
    mode = body[2 * SIZEOF_NETWORK_LONG + 1];
    flags = body[2 * SIZEOF_NETWORK_LONG + 2];
    const unsigned char *verticesBody = body + fixedBytes;
    const unsigned char *facesBody =
	verticesBody + vertices * vertexStride;
    compact_coverage_vertex_bounds vertexBounds;
    const bool parallelBounds =
	vertices >= compact_coverage_parallel_vertex_threshold &&
	compact_coverage_try_acquire_parallel_vertex_scan();
    if (parallelBounds) {
	const size_t workerCount = std::min(
	    compact_coverage_parallel_vertex_workers,
	    std::max<size_t>(1, bu_avail_cpus()));
	std::array<compact_coverage_vertex_bounds,
	    compact_coverage_parallel_vertex_workers> partialBounds;
	std::array<std::thread, compact_coverage_parallel_vertex_workers> workers;
	size_t launched = 0;
	try {
	    for (size_t worker = 0; worker < workerCount; ++worker) {
		const size_t first = vertices / workerCount * worker;
		const size_t last = vertices / workerCount * (worker + 1);
		workers[worker] = std::thread([&, worker, first, last]() {
		    partialBounds[worker] = compact_coverage_scan_serialized_vertices(
			verticesBody, first, last - first, vertexStride);
		});
		launched++;
	    }
	} catch (const std::exception &) {
	}
	for (size_t worker = 0; worker < launched; ++worker)
	    workers[worker].join();
	compact_coverage_active_vertex_scans.fetch_sub(1);
	if (launched == workerCount) {
	    for (size_t worker = 0; worker < workerCount; ++worker)
		compact_coverage_merge_vertex_bounds(vertexBounds,
		    partialBounds[worker]);
	} else {
	    vertexBounds = compact_coverage_scan_serialized_vertices(
		verticesBody, 0, vertices, vertexStride);
	}
    } else {
	vertexBounds = compact_coverage_scan_serialized_vertices(
	    verticesBody, 0, vertices, vertexStride);
    }
    if (vertexBounds.finiteVertices && faces) {
	const size_t sampleCount = std::min<size_t>(32, faces);
	uint64_t hash = 1469598103934665603ULL;
	const auto mix = [&hash](uint64_t word) {
	    hash ^= word;
	    hash *= 1099511628211ULL;
	};
	mix(static_cast<uint64_t>(vertices));
	mix(static_cast<uint64_t>(faces));
	bool valid = true;
	for (size_t sample = 0; sample < sampleCount && valid; ++sample) {
	    const size_t face = sampleCount > 1 ?
		sample * (faces - 1) / (sampleCount - 1) : 0;
	    const unsigned char *indices = facesBody + face * faceStride;
	    uint32_t vertex[3] = {
		static_cast<uint32_t>(BU_GLONG(indices)),
		static_cast<uint32_t>(BU_GLONG(
		    indices + SIZEOF_NETWORK_LONG)),
		static_cast<uint32_t>(BU_GLONG(
		    indices + 2 * SIZEOF_NETWORK_LONG))
	    };
	    for (size_t edge = 0; edge < 3; ++edge) {
		const size_t first = vertex[edge];
		const size_t second = vertex[(edge + 1) % 3];
		if (first >= vertices || second >= vertices) {
		    valid = false;
		    break;
		}
		double a[ELEMENTS_PER_POINT];
		double b[ELEMENTS_PER_POINT];
		bu_cv_ntohd(reinterpret_cast<unsigned char *>(a),
		    verticesBody + first * vertexStride,
		    ELEMENTS_PER_POINT);
		bu_cv_ntohd(reinterpret_cast<unsigned char *>(b),
		    verticesBody + second * vertexStride,
		    ELEMENTS_PER_POINT);
		const double dx = a[X] - b[X];
		const double dy = a[Y] - b[Y];
		const double dz = a[Z] - b[Z];
		const double lengthSquared = dx * dx + dy * dy + dz * dz;
		if (!std::isfinite(lengthSquared) || lengthSquared < 0.0) {
		    valid = false;
		    break;
		}
		int64_t quantized = std::numeric_limits<int64_t>::min();
		if (lengthSquared > SMALL_FASTF) {
		    const double scaled = log2(lengthSquared) * 262144.0;
		    if (!std::isfinite(scaled) ||
			scaled > static_cast<double>(
			    std::numeric_limits<int64_t>::max()) ||
			scaled < static_cast<double>(
			    std::numeric_limits<int64_t>::min())) {
			valid = false;
			break;
		    }
		    quantized = static_cast<int64_t>(llround(scaled));
		}
		mix(static_cast<uint64_t>(quantized));
	    }
	}
	if (valid) {
	    sampleFingerprint = hash;
	    sampleFingerprintValid = true;
	}
    }
    if (!borrowed)
	bu_free_external(&external);
    if (!vertexBounds.finiteVertices)
	return CompactCoverageBotStatus::INVALID;

    bounds = SbBox3f(
	SbVec3f(static_cast<float>(vertexBounds.minimum[X]),
	    static_cast<float>(vertexBounds.minimum[Y]),
	    static_cast<float>(vertexBounds.minimum[Z])),
	SbVec3f(static_cast<float>(vertexBounds.maximum[X]),
	    static_cast<float>(vertexBounds.maximum[Y]),
	    static_cast<float>(vertexBounds.maximum[Z])));
    return bounds.isEmpty() ? CompactCoverageBotStatus::INVALID :
	CompactCoverageBotStatus::READY;
}

struct compact_coverage_mapped_database {
    struct bu_mapped_file *file = NULL;
    ~compact_coverage_mapped_database(void)
    {
	if (file)
	    bu_close_mapped_file(file);
    }
};

/* The mapped coverage path is a latency optimization for ordinary databases,
 * not a residency contract.  Mapping a multi-gigabyte editable database a
 * second time during startup competes with the directory, worker, cache, and
 * renderer working sets.  Large sources already have a bounded per-object
 * external-view fallback, which is the correct path when whole-file virtual
 * reservation would make an otherwise cheap coverage pass fail. */
static bool
compact_coverage_can_map_database(const struct db_i *dbip)
{
    static const int maximumMappedDatabaseBytes = 512 * 1024 * 1024;
    if (!dbip || !dbip->dbi_filename || !dbip->dbi_filename[0])
	return false;
    const int databaseBytes = bu_file_size(dbip->dbi_filename);
    return databaseBytes >= 0 &&
	databaseBytes <= maximumMappedDatabaseBytes;
}

static bool
compact_coverage_serialized_point(
	const BObolSerializedBotView &view, size_t index,
	point_t point)
{
    if (!view.vertices || index >= view.vertexCount)
	return false;
    double decoded[ELEMENTS_PER_POINT];
    bu_cv_ntohd(reinterpret_cast<unsigned char *>(decoded),
	view.vertices + index * view.vertexStride, ELEMENTS_PER_POINT);
    if (!std::isfinite(decoded[X]) || !std::isfinite(decoded[Y]) ||
	!std::isfinite(decoded[Z]))
	return false;
    VMOVE(point, decoded);
    return true;
}

/* Prove that candidate is a rigid transform of representative without
 * allocating either authored mesh.  Xpush preserves vertex and face order,
 * so corresponding well-separated vertices define the candidate transform;
 * every vertex and the complete topology must then verify.  A reordered or
 * numerically ambiguous mesh is a safe false negative and loads independently.
 */
static bool
compact_coverage_serialized_rigid_match(struct db_i *dbip,
	struct directory *representativeDp, struct directory *candidateDp,
	SbMatrix &representativeToCandidate,
	const struct bu_mapped_file *mappedFile)
{
    const bool verbose = getenv("BOBOL_DRAW_TIMING_VERBOSE") != NULL;
    const char *representativeName = representativeDp &&
	representativeDp->d_namep ? representativeDp->d_namep : "?";
    const char *candidateName = candidateDp && candidateDp->d_namep ?
	candidateDp->d_namep : "?";
#define COVERAGE_REUSE_REJECT(reason) do { \
    if (verbose) \
        bu_log("[obol-timing] serialized reuse reject: %s -> %s: %s\n", \
            candidateName, representativeName, reason); \
    return false; \
} while (0)
    representativeToCandidate = SbMatrix::identity();
    BObolSerializedBotView representative;
    BObolSerializedBotView candidate;
    if (!bobol_serialized_bot_view(dbip, representativeDp,
	representative, mappedFile) ||
	!bobol_serialized_bot_view(dbip, candidateDp, candidate, mappedFile))
	COVERAGE_REUSE_REJECT("serialized view");
    if (representative.vertexCount != candidate.vertexCount ||
	representative.faceCount != candidate.faceCount)
	COVERAGE_REUSE_REJECT("count mismatch");
    if (memcmp(representative.faces, candidate.faces,
	representative.faceCount * representative.faceStride) != 0)
	COVERAGE_REUSE_REJECT("topology mismatch");

    point_t rp0;
    point_t cp0;
    if (!compact_coverage_serialized_point(representative, 0, rp0) ||
	!compact_coverage_serialized_point(candidate, 0, cp0))
	COVERAGE_REUSE_REJECT("invalid first point");

    size_t p1Index = 0;
    double farthestSquared = 0.0;
    for (size_t i = 1; i < representative.vertexCount; ++i) {
	point_t point;
	if (!compact_coverage_serialized_point(representative, i, point))
	    COVERAGE_REUSE_REJECT("invalid representative point");
	vect_t delta;
	VSUB2(delta, point, rp0);
	const double lengthSquared = MAGSQ(delta);
	if (lengthSquared > farthestSquared) {
	    farthestSquared = lengthSquared;
	    p1Index = i;
	}
    }
    if (p1Index == 0 || farthestSquared <= SMALL_FASTF)
	COVERAGE_REUSE_REJECT("degenerate primary axis");

    point_t rp1;
    point_t cp1;
    if (!compact_coverage_serialized_point(representative, p1Index, rp1) ||
	!compact_coverage_serialized_point(candidate, p1Index, cp1))
	COVERAGE_REUSE_REJECT("invalid primary correspondence");
    vect_t rx;
    vect_t cx;
    VSUB2(rx, rp1, rp0);
    VSUB2(cx, cp1, cp0);
    const double characteristicExtent = MAGNITUDE(rx);
    const double candidateExtent = MAGNITUDE(cx);
    const double tolerance = std::max(static_cast<double>(VUNITIZE_TOL),
	characteristicExtent * 1.0e-9);
    if (fabs(characteristicExtent - candidateExtent) > tolerance)
	COVERAGE_REUSE_REJECT("primary extent mismatch");
    VUNITIZE(rx);
    VUNITIZE(cx);

    size_t p2Index = 0;
    double farthestNormalSquared = 0.0;
    for (size_t i = 1; i < representative.vertexCount; ++i) {
	point_t point;
	if (!compact_coverage_serialized_point(representative, i, point))
	    COVERAGE_REUSE_REJECT("invalid normal search point");
	vect_t delta;
	vect_t normal;
	VSUB2(delta, point, rp0);
	VCROSS(normal, rx, delta);
	const double lengthSquared = MAGSQ(normal);
	if (lengthSquared > farthestNormalSquared) {
	    farthestNormalSquared = lengthSquared;
	    p2Index = i;
	}
    }
    if (p2Index == 0 ||
	farthestNormalSquared <= tolerance * tolerance)
	COVERAGE_REUSE_REJECT("degenerate secondary axis");

    point_t rp2;
    point_t cp2;
    if (!compact_coverage_serialized_point(representative, p2Index, rp2) ||
	!compact_coverage_serialized_point(candidate, p2Index, cp2))
	COVERAGE_REUSE_REJECT("invalid secondary correspondence");
    vect_t rdelta;
    vect_t cdelta;
    vect_t rz;
    vect_t cz;
    vect_t ry;
    vect_t cy;
    VSUB2(rdelta, rp2, rp0);
    VSUB2(cdelta, cp2, cp0);
    VCROSS(rz, rx, rdelta);
    VCROSS(cz, cx, cdelta);
    if (MAGNITUDE(rz) <= tolerance || MAGNITUDE(cz) <= tolerance)
	COVERAGE_REUSE_REJECT("secondary extent mismatch");
    VUNITIZE(rz);
    VUNITIZE(cz);
    VCROSS(ry, rz, rx);
    VCROSS(cy, cz, cx);
    VUNITIZE(ry);
    VUNITIZE(cy);

    struct bg_pca_frame sourceFrame = {};
    struct bg_pca_frame targetFrame = {};
    VMOVE(sourceFrame.center, rp0);
    VMOVE(sourceFrame.xaxis, rx);
    VMOVE(sourceFrame.yaxis, ry);
    VMOVE(sourceFrame.zaxis, rz);
    VMOVE(targetFrame.center, cp0);
    VMOVE(targetFrame.xaxis, cx);
    VMOVE(targetFrame.yaxis, cy);
    VMOVE(targetFrame.zaxis, cz);
    VSETALL(sourceFrame.singular_values, 1.0);
    VSETALL(targetFrame.singular_values, 1.0);
    mat_t relative;
    if (bg_pca_frame_relative_matrix(relative, &sourceFrame,
	&targetFrame) != BRLCAD_OK)
	COVERAGE_REUSE_REJECT("relative matrix");

    static const size_t vertexBlock = 256;
    double sourcePoints[vertexBlock * ELEMENTS_PER_POINT];
    double targetPoints[vertexBlock * ELEMENTS_PER_POINT];
    for (size_t first = 0; first < representative.vertexCount;
	 first += vertexBlock) {
	const size_t count = std::min(vertexBlock,
	    representative.vertexCount - first);
	bu_cv_ntohd(reinterpret_cast<unsigned char *>(sourcePoints),
	    representative.vertices + first * representative.vertexStride,
	    count * ELEMENTS_PER_POINT);
	bu_cv_ntohd(reinterpret_cast<unsigned char *>(targetPoints),
	    candidate.vertices + first * candidate.vertexStride,
	    count * ELEMENTS_PER_POINT);
	for (size_t i = 0; i < count; ++i) {
	    point_t transformed;
	    MAT4X3PNT(transformed, relative,
		sourcePoints + i * ELEMENTS_PER_POINT);
	    const double *target = targetPoints + i * ELEMENTS_PER_POINT;
	    if (!std::isfinite(target[X]) || !std::isfinite(target[Y]) ||
		!std::isfinite(target[Z]) ||
		fabs(transformed[X] - target[X]) > tolerance ||
		fabs(transformed[Y] - target[Y]) > tolerance ||
		fabs(transformed[Z] - target[Z]) > tolerance)
		COVERAGE_REUSE_REJECT("all-vertex verification");
	}
    }
    representativeToCandidate = mat_to_sbmatrix(relative);
#undef COVERAGE_REUSE_REJECT
    return true;
}

static bool
compact_coverage_primitive_bounds(struct db_i *dbip, struct directory *dp,
	SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!dbip || !dp)
	return false;
    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, dbip, NULL) < 0)
	return false;
    const bool valid = local_bounds_from_internal(&intern, bounds);
    rt_db_free_internal(&intern);
    return valid && !bounds.isEmpty();
}

static size_t
compact_coverage_bot_source_bytes(const struct rt_bot_internal *bot)
{
    if (!bot)
	return 0;
    size_t bytes = sizeof(*bot);
    const auto add = [&bytes](size_t count, size_t elementSize) {
	if (bytes == SIZE_MAX || !count || !elementSize)
	    return;
	if (count > (SIZE_MAX - bytes) / elementSize)
	    bytes = SIZE_MAX;
	else
	    bytes += count * elementSize;
    };
    add(bot->num_vertices, 3 * sizeof(fastf_t));
    add(bot->num_faces, 3 * sizeof(int));
    if (bot->thickness)
	add(bot->num_faces, sizeof(fastf_t));
    add(bot->num_normals, 3 * sizeof(fastf_t));
    add(bot->num_face_normals, 3 * sizeof(int));
    add(bot->num_uvs, 3 * sizeof(fastf_t));
    add(bot->num_face_uvs, 3 * sizeof(int));
    return bytes;
}

static union tree *
compact_coverage_collect_leaf(struct db_tree_state *tsp,
	const struct db_full_path *pathp, struct directory *dp,
	void *clientData)
{
    compact_coverage_collect *collect =
	static_cast<compact_coverage_collect *>(clientData);
    if (!collect || !collect->source || !collect->cache ||
	!collect->stream || !tsp || !tsp->ts_dbip || !pathp || !dp)
	return TREE_NULL;
    if (collect->stream->isCancelled())
	return make_nop_tree();
    const std::string instanceIdentity =
	realize_walk_instance_identity(tsp, pathp);
    if (instanceIdentity.empty() ||
	!collect->seenInstances.insert(instanceIdentity).second)
	return make_nop_tree();
    const std::string occurrenceIdentity =
	realize_walk_occurrence_identity(pathp);
    const uint32_t duplicateOrdinal =
	collect->occurrenceCounts[occurrenceIdentity]++;

    char *rawPath = db_path_to_string(pathp);
    if (!rawPath || !rawPath[0]) {
	if (rawPath)
	    bu_free(rawPath, "compact coverage path");
	return make_nop_tree();
    }

    SbBox3f unusedBounds;
    std::string cacheKey = realize_geometry_cache_key(dp);
    source_lod_cache_key_append(cacheKey, collect->source, unusedBounds,
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP);
    size_t assetIndex = 0;
    auto foundAsset = collect->assetIndices.find(dp);
    if (foundAsset == collect->assetIndices.end()) {
	assetIndex = collect->assets.size();
	std::unique_ptr<compact_coverage_asset> asset(
	    new compact_coverage_asset);
	asset->dp = dp;
	asset->cacheKey = cacheKey;
	asset->assetPath = rawPath;
	asset->estimatedWorkingSetBytes =
	    compact_coverage_working_set_estimate(tsp->ts_dbip, dp);
	collect->assets.push_back(std::move(asset));
	collect->assetIndices[dp] = assetIndex;
    } else {
	assetIndex = foundAsset->second;
    }

    compact_coverage_occurrence occurrence;
    occurrence.localTransform = mat_to_sbmatrix(tsp->ts_mat);
    occurrence.assetPath = rawPath;
    std::string semanticPath = rawPath;
    if (duplicateOrdinal > 0) {
	char suffix[32] = {0};
	snprintf(suffix, sizeof(suffix), "@%u", duplicateOrdinal);
	semanticPath += suffix;
    }
    const char *sourceType =
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP ? "brep" :
	(dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT ? "bot" :
	 "primitive");
    occurrence.summary = compact_occurrence_tree_summary(
	collect->source, tsp, pathp, semanticPath.c_str(), dp->d_namep,
	sourceType,
	"aabb", collect->revision,
	BObolRealizedShapeSummary::SHAPE_MESH, collect->materialSweep);
    occurrence.occurrenceIndex = pathp->fp_cinst && pathp->fp_len ?
	static_cast<uint32_t>(DB_FULL_PATH_GET_COMB_INST(pathp,
		pathp->fp_len - 1)) : 0;
    occurrence.booleanOperation = (tsp->ts_sofar & TS_SOFAR_MINUS) ?
	SoBRLDatabaseSource::BOOLEAN_SUBTRACT :
	((tsp->ts_sofar & TS_SOFAR_INTER) ?
	 SoBRLDatabaseSource::BOOLEAN_INTERSECT :
	 SoBRLDatabaseSource::BOOLEAN_UNION);
    if (occurrence.booleanOperation !=
	    SoBRLDatabaseSource::BOOLEAN_UNION ||
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP)
	collect->aggregateBoundsExact = false;
    compact_coverage_work_item work;
    work.asset = collect->assets[assetIndex].get();
    work.occurrence = std::move(occurrence);
    collect->producerWork.push_back(std::move(work));
    collect->occurrenceCount++;
    /* The first leaf is latency-critical: it can supply both the expanding
     * scene overview and the first useful leaf box while the hierarchy walk
     * continues.  Waiting for a full producer batch accidentally made a
     * one-leaf database (and any slow first branch with fewer than 16 leaves)
     * wait for the complete walk before a coverage worker could start.  Keep
     * batching for throughput after that first publication. */
    if (collect->occurrenceCount == 1 ||
	collect->producerWork.size() >= compact_coverage_work_batch_size) {
	{
	    std::lock_guard<std::mutex> guard(collect->workMutex);
	    collect->work.insert(collect->work.end(),
		std::make_move_iterator(collect->producerWork.begin()),
		std::make_move_iterator(collect->producerWork.end()));
	    collect->producerWork.clear();
	}
	collect->workReady.notify_one();
    }
    bu_free(rawPath, "compact coverage path");
    return make_nop_tree();
}

static std::string
compact_coverage_broad_key(const compact_coverage_asset &asset)
{
    if (!asset.coverageReady || !asset.sampleFingerprintValid ||
	asset.flags != 0 || asset.mode == RT_BOT_PLATE ||
	asset.mode == RT_BOT_PLATE_NOCOS)
	return std::string();
    char key[192] = {0};
	snprintf(key, sizeof(key), "%zu:%zu:%u:%u:%u:%016llx",
	asset.vertexCount, asset.faceCount,
	static_cast<unsigned int>(asset.mode),
	static_cast<unsigned int>(asset.orientation),
	static_cast<unsigned int>(asset.flags),
	static_cast<unsigned long long>(asset.sampleFingerprint));
    return std::string(key);
}

static bool
compact_coverage_lod_asset_oriented_bounds(
	struct db_i *dbip, const BObolDrawLodAssetRecord &mapping,
	std::array<SbVec3f, 8> &objectBounds)
{
    if (!dbip || !mapping.assetName[0])
	return false;

    const BObolDrawLodAssetRecord *metadata = &mapping;
    BObolDrawLodAssetRecord canonical;
    if (mapping.assetOrientedBoundsValid != 1) {
	if (bobol_draw_lod_asset_cache_get(dbip, mapping.assetName,
		&canonical) != BRLCAD_OK ||
	    canonical.assetOrientedBoundsValid != 1 ||
	    !BU_STR_EQUAL(canonical.assetName, mapping.assetName))
	    return false;
	metadata = &canonical;
    }

    const SbMatrix assetToObject = mat_to_sbmatrix(mapping.assetToObject);
    for (size_t corner = 0; corner < objectBounds.size(); ++corner) {
	const point_t &point = metadata->assetOrientedBounds[corner];
	const SbVec3f source(
	    static_cast<float>(point[X]), static_cast<float>(point[Y]),
	    static_cast<float>(point[Z]));
	assetToObject.multVecMatrix(source, objectBounds[corner]);
	if (!std::isfinite(objectBounds[corner][0]) ||
	    !std::isfinite(objectBounds[corner][1]) ||
	    !std::isfinite(objectBounds[corner][2]))
	    return false;
    }
    return true;
}

static std::shared_ptr<const Obol::PartGeometry>
compact_coverage_overview_geometry(const SbBox3f &bounds)
{
    Obol::PartGeometryBuilder geometry;
    if (!cad_wire_part_geometry_from_aabb(bounds, geometry))
	return std::shared_ptr<const Obol::PartGeometry>();
    /* The overview is already one aggregate for the complete target.  It must
     * remain a visible extent, not enter the leaf subpixel-point classifier. */
    geometry.subpixelProxyEligible = false;
    geometry.structuralProxy = true;
    return bobol_cad_build_geometry(
	std::move(geometry), "coverage overview");
}

static BObolCompactOccurrence
compact_coverage_leaf_occurrence(const compact_coverage_asset &asset,
	const compact_coverage_occurrence &leaf, bool includeSourceRequest)
{
    BObolCompactOccurrence occurrence;
    occurrence.geometry = asset.coverageGeometry;
    occurrence.geometryTransform = asset.proxyGeometryTransform;
    occurrence.localTransform = leaf.localTransform;
    occurrence.lodBacked = asset.lodEligible ? TRUE : FALSE;
    occurrence.summary = leaf.summary;
    occurrence.sourceMeshRequestValid =
	includeSourceRequest && asset.sourceMeshReady ? TRUE : FALSE;
    if (occurrence.sourceMeshRequestValid) {
	occurrence.sourceMeshRequest = asset.sourceMeshRequest;
	compact_source_mesh_request_sync(
	    occurrence.sourceMeshRequest, occurrence.summary);
	compact_summary_lod_from_source_mesh_request(
	    occurrence.summary, occurrence.sourceMeshRequest);
    }
    occurrence.occurrenceIndex = leaf.occurrenceIndex;
    occurrence.booleanOperation = leaf.booleanOperation;
    occurrence.summary.lodAvailable = asset.lodEligible ? TRUE : FALSE;
    occurrence.summary.lodActiveCut = BOBOL_LOD_QUALITY_PROXY;
    occurrence.summary.lodFaceCount = asset.faceCount;
    occurrence.summary.lodPointCount = asset.vertexCount;
    occurrence.summary.lodOriginalPointCount = asset.vertexCount;
    occurrence.summary.lodBoundsMin = asset.coverageBounds.getMin();
    occurrence.summary.lodBoundsMax = asset.coverageBounds.getMax();
    /* The source request is a valid lazy refinement contract, but its current
     * presentation is an AABB rather than a resident PoP cut.  Do not encode
     * proxy quality enums and authored source counts in the provider's
     * resident-prefix fields: the submitter would treat them as a warm full
     * mesh and reserve whole-source conversion scratch before the cold
     * spatial producer can run. */
    if (occurrence.sourceMeshRequestValid) {
	occurrence.sourceMeshRequest.lodAvailable = 0;
	occurrence.sourceMeshRequest.lodActiveCut = -1;
	occurrence.sourceMeshRequest.lodFaceCount = 0;
	occurrence.sourceMeshRequest.lodPointCount = 0;
	occurrence.sourceMeshRequest.lodOriginalPointCount = 0;
	occurrence.sourceMeshRequest.lodBoundsMin =
	    asset.coverageBounds.getMin();
	occurrence.sourceMeshRequest.lodBoundsMax =
	    asset.coverageBounds.getMax();
    }
    occurrence.summary.boundsValid = TRUE;
    occurrence.summary.bounds = asset.coverageBounds;
    return occurrence;
}

/* Rehydrate only the terminal subset of a complete warm manifest.  Lazy mesh
 * records are already sufficient for view-driven LoD and must not be scanned
 * again.  Analytic/BREP records retain their semantic occurrence state in the
 * manifest, so this seed recreates the ordinary detail-work inputs without a
 * second hierarchy walk. */
static bool
compact_coverage_seed_warm_terminal_work(
	SoBRLDatabaseSource *source,
	BObolCompactOccurrenceStream *stream,
	compact_coverage_collect &collect,
	uint32_t revision)
{
    if (!source || !stream || !source->getDatabase() ||
	!stream->hasWarmCensusComplete() ||
	stream->hasWarmCoverageComplete())
	return false;

    std::vector<BObolCompactManifestOccurrence> records;
    if (!stream->takeWarmTerminalOccurrences(records) || records.empty())
	return false;
    const auto rejectWarmCensus = [stream]() {
	stream->setWarmCensusComplete(false);
	return false;
    };

    /* Validate the complete replay set before publishing it to collect.  A
     * stale or malformed record falls back to the ordinary hierarchy walk;
     * it must not leave partial assets or duplicate detail work behind. */
    std::vector<std::unique_ptr<compact_coverage_asset>> stagedAssets;
    stagedAssets.reserve(records.size());
    std::deque<compact_coverage_work_item> stagedWork;
    std::unordered_map<struct directory *, compact_coverage_asset *>
	assetsByDirectory;
    assetsByDirectory.reserve(records.size());
    for (const BObolCompactManifestOccurrence &record : records) {
	const char *sourceName = record.sourceName.getString();
	const char *path = record.path.getString();
	if (!sourceName || !sourceName[0] || !path || !path[0] ||
	    record.bounds.isEmpty() || record.sourceMeshRequestValid)
	    return rejectWarmCensus();
	struct directory *dp = db_lookup(source->getDatabase(), sourceName,
	    LOOKUP_QUIET);
	if (dp == RT_DIR_NULL ||
	    dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT)
	    return rejectWarmCensus();

	compact_coverage_asset *asset = NULL;
	auto found = assetsByDirectory.find(dp);
	if (found == assetsByDirectory.end()) {
	    std::unique_ptr<compact_coverage_asset> owned(
		new compact_coverage_asset);
	    owned->dp = dp;
	    owned->assetPath = sourceName;
	    owned->estimatedWorkingSetBytes =
		compact_coverage_working_set_estimate(
		    source->getDatabase(), dp);
	    owned->coverageBounds = record.bounds;
	    owned->coverageGeometry = bobol_cad_structural_bounds_geometry(
		record.bounds, owned->proxyGeometryTransform);
	    owned->coverageReady = owned->coverageGeometry ? true : false;
	    asset = owned.get();
	    stagedAssets.push_back(std::move(owned));
	    assetsByDirectory.emplace(dp, asset);
	} else {
	    asset = found->second;
	}
	if (!asset || !asset->coverageReady)
	    return rejectWarmCensus();

	compact_coverage_occurrence occurrence;
	occurrence.localTransform = record.localTransform;
	occurrence.assetPath = sourceName;
	occurrence.occurrenceIndex = record.occurrenceIndex;
	occurrence.booleanOperation = record.booleanOperation;
	occurrence.summary = compact_occurrence_summary(source, path,
	    sourceName, "primitive", "aabb", revision,
	    BObolRealizedShapeSummary::SHAPE_MESH);
	occurrence.summary.regionId = record.regionId;
	occurrence.summary.airCode = record.airCode;
	occurrence.summary.materialId = record.materialId;
	occurrence.summary.los = record.los;
	occurrence.summary.materialColorValid = record.materialColorValid;
	occurrence.summary.materialColor = record.materialColor;
	occurrence.summary.materialShader = record.materialShader;
	occurrence.summary.boundsValid = TRUE;
	occurrence.summary.bounds = record.bounds;

	compact_coverage_work_item item;
	item.asset = asset;
	item.occurrence = std::move(occurrence);
	stagedWork.push_back(std::move(item));
    }

    const size_t occurrenceCount = stream->getExpectedCount();
    if (occurrenceCount == 0 || stagedWork.empty())
	return rejectWarmCensus();
    collect.assets = std::move(stagedAssets);
    collect.detailWork = std::move(stagedWork);
    collect.occurrenceCount = occurrenceCount;
    return true;
}

static int
compact_stream_publish_parallel_coverage(
	SoBRLDatabaseSource *source,
	BObolDatabaseSourceRealizationCache *cache,
	const char *treeName,
	BObolCompactOccurrenceStream *stream,
	uint32_t revision,
	bool *authoritativeStreamOut)
{
    if (authoritativeStreamOut)
	*authoritativeStreamOut = false;
    if (!source || !cache || !treeName || !treeName[0] || !stream ||
	!source->getDatabase() || source->lodBotThreshold.getValue() == 0)
	return 0;

    const int64_t collectStart = bu_gettime();
    BObolMaterialColorSweep materialSweep(source->getDatabase());
    compact_coverage_collect collect;
    collect.source = source;
    collect.cache = cache;
    collect.stream = stream;
    collect.materialSweep = &materialSweep;
    collect.revision = revision;
    const bool selectiveWarmReplay =
	compact_coverage_seed_warm_terminal_work(source, stream, collect,
	    revision);
    compact_coverage_mapped_database mappedDatabase;
    if (!selectiveWarmReplay &&
	compact_coverage_can_map_database(source->getDatabase())) {
	mappedDatabase.file = bu_open_mapped_file(
	    source->getDatabase()->dbi_filename,
	    "obol-coverage-reuse-v1");
    }
    std::atomic<size_t> publishedBoxes(0);
    std::atomic<size_t> publishedTerminalOnly(0);
    std::atomic<size_t> publishedContracts(0);
    std::atomic<size_t> terminalEmptyOccurrences(0);
    std::mutex aggregateMutex;
    SbBox3f aggregateBounds;
    aggregateBounds.makeEmpty();
    SbBox3f certifiedWarmBounds;
    const bool haveCertifiedWarmBounds =
	stream->hasCoverageBoundsComplete() &&
	stream->getCoverageBounds(certifiedWarmBounds) &&
	!certifiedWarmBounds.isEmpty();
    if (selectiveWarmReplay)
	publishedBoxes.store(collect.occurrenceCount);
    int64_t lastOverviewPublication = 0;
    size_t workerCount = bu_avail_cpus();
    workerCount = std::max<size_t>(1, std::min<size_t>(workerCount, 32));
    collect.producerWork.reserve(compact_coverage_work_batch_size);
    const auto coverageWorker = [&]() {
	/* Bounds and the first useful overview are foreground draw latency, not
	 * speculative refinement.  Lowering these workers' priority let unrelated
	 * renderer/cache work starve a cold scene into several blank seconds.
	 * Full-array import below remains niced. */
	std::vector<compact_coverage_work_item> items;
	items.reserve(compact_coverage_work_batch_size);
	std::vector<BObolCompactOccurrence> publications;
	publications.reserve(compact_coverage_work_batch_size + 1);
	for (;;) {
	    items.clear();
	    publications.clear();
	    {
		std::unique_lock<std::mutex> lock(collect.workMutex);
		collect.workReady.wait(lock, [&]() {
		    return collect.producerDone || !collect.work.empty() ||
			stream->isCancelled();
		});
		if (stream->isCancelled())
		    break;
		if (collect.work.empty()) {
		    if (collect.producerDone)
			break;
		    continue;
		}
		const size_t count = std::min(compact_coverage_work_batch_size,
		    collect.work.size());
		for (size_t i = 0; i < count; ++i) {
		    items.push_back(std::move(collect.work.front()));
		    collect.work.pop_front();
		}
	    }
	    for (compact_coverage_work_item &item : items) {
	    if (!item.asset)
		continue;
	    compact_coverage_asset &asset = *item.asset;
	    std::call_once(asset.coverageOnce, [&]() {
		size_t coverageWorkingSet = asset.estimatedWorkingSetBytes;
		if (db_version(source->getDatabase()) == 5 && asset.dp &&
		    asset.dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) {
		    const size_t fixedBytes = 1024ULL * 1024ULL;
		    size_t borrowedBytes = 0;
		    const bool zeroCopy = mappedDatabase.file || db_external_view(
			source->getDatabase(), asset.dp,
			&borrowedBytes) != NULL;
		    const size_t encodedBytes = zeroCopy ? 0 : asset.dp->d_len;
		    coverageWorkingSet =
			encodedBytes > SIZE_MAX - fixedBytes ?
			SIZE_MAX : encodedBytes + fixedBytes;
		}
		if (!bobol_lod_working_set_acquire(coverageWorkingSet))
		    return;
		if (!stream->isCancelled() && asset.dp &&
		    asset.dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) {
		    const CompactCoverageBotStatus coverageStatus =
			compact_coverage_bot_bounds(
			source->getDatabase(), asset.dp,
			asset.coverageBounds, asset.vertexCount,
			asset.faceCount, asset.mode, asset.orientation,
			asset.flags, asset.sampleFingerprint,
			asset.sampleFingerprintValid,
			mappedDatabase.file);
		    asset.coverageReady =
			coverageStatus == CompactCoverageBotStatus::READY;
		    asset.terminalEmpty =
			coverageStatus == CompactCoverageBotStatus::EMPTY;
		    asset.lodEligible = asset.coverageReady &&
			asset.faceCount >=
			    source->lodBotThreshold.getValue();
		} else if (!stream->isCancelled() && asset.dp) {
		    asset.coverageReady = compact_coverage_primitive_bounds(
			source->getDatabase(), asset.dp,
			asset.coverageBounds);
		    asset.lodEligible = asset.coverageReady &&
			asset.dp->d_minor_type ==
			    DB5_MINORTYPE_BRLCAD_BREP;
		}
		if (asset.coverageReady) {
		    asset.coverageGeometry = bobol_cad_structural_bounds_geometry(
			asset.coverageBounds, asset.proxyGeometryTransform);
		    asset.coverageReady =
			asset.coverageGeometry ? true : false;
		    /*
		     * A BoT's compact leaf contract needs identity, tight bounds and
		     * counts; it does not need the full authored arrays.  Publish that
		     * contract with the standing leaf box as soon as the v5 coverage
		     * decoder has supplied those facts.  The view-LoD provider can then
		     * open a warm PoP payload, or import/generate a cold one, only after
		     * projection and the scene allocator have admitted this occurrence.
		     *
		     * The previous second phase unconditionally decoded every BoT before
		     * any view decision.  On a 150k-part scene that made subpixel and
		     * offscreen leaves consume the same CPU, cache I/O and resident-memory
		     * path as prominent visible parts, preventing convergence even with a
		     * theoretically warm cache.
		     */
		    if (asset.coverageReady && asset.lodEligible && asset.dp &&
			asset.dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) {
			asset.sourceMeshRequest.clear();
			asset.sourceMeshRequest.meshAssetPath =
			    asset.assetPath.c_str();
			asset.sourceMeshRequest.meshAssetName =
			    asset.dp->d_namep ? asset.dp->d_namep : "";
			asset.sourceMeshRequest.meshAssetBounds =
			    asset.coverageBounds;
			asset.sourceMeshRequest.bounds = asset.coverageBounds;
			asset.sourceMeshRequest.faceCount = asset.faceCount;
			asset.sourceMeshRequest.pointCount = asset.vertexCount;
			asset.sourceMeshRequest.meshAssetTransform =
			    SbMatrix::identity();
			asset.sourceMeshReady = true;
			/* Exact transformed mappings are validated against both database
			 * objects by the draw cache.  Reapply one here so a warm draw can
			 * publish the canonical contract immediately and skip both the
			 * all-vertex proof and a duplicate PoP request. */
			BObolDrawLodAssetRecord cachedAsset;
			if (bobol_draw_lod_asset_cache_get(
				source->getDatabase(), asset.dp->d_namep,
				&cachedAsset) == BRLCAD_OK &&
			    cachedAsset.faceCount == asset.faceCount &&
			    cachedAsset.pointCount == asset.vertexCount) {
			    struct directory *assetDp = db_lookup(
				source->getDatabase(), cachedAsset.assetName,
				LOOKUP_QUIET);
			    if (assetDp != RT_DIR_NULL) {
				asset.sourceMeshRequest.meshAssetPath =
				    cachedAsset.assetName;
				asset.sourceMeshRequest.meshAssetName =
				    cachedAsset.assetName;
				asset.sourceMeshRequest.meshAssetBounds =
				    SbBox3f(SbVec3f(
					cachedAsset.assetBoundsMin[X],
					cachedAsset.assetBoundsMin[Y],
					cachedAsset.assetBoundsMin[Z]),
					SbVec3f(cachedAsset.assetBoundsMax[X],
					    cachedAsset.assetBoundsMax[Y],
					    cachedAsset.assetBoundsMax[Z]));
				asset.sourceMeshRequest.meshAssetTransform =
				    mat_to_sbmatrix(cachedAsset.assetToObject);
				asset.sourceMeshMappingCached = true;
				asset.deferSourceMeshContract = false;
				std::array<SbVec3f, 8> orientedBounds;
				if (compact_coverage_lod_asset_oriented_bounds(
					source->getDatabase(), cachedAsset,
					orientedBounds)) {
				    std::shared_ptr<const Obol::PartGeometry>
					orientedCoverage =
					    bobol_cad_structural_bounds_geometry(
						asset.coverageBounds,
						asset.proxyGeometryTransform,
						orientedBounds.data());
				    if (orientedCoverage)
					asset.coverageGeometry =
					    std::move(orientedCoverage);
				}
			    }
			}
			/* The first member of a rigid-invariant group may begin loading
			 * immediately.  Later members publish their boxes now but defer
			 * their source contracts until the exact serialized proof below,
			 * preventing duplicate giant PoP builds from racing that proof. */
			const std::string broadKey =
			    compact_coverage_broad_key(asset);
			if (!broadKey.empty() &&
			    !asset.sourceMeshMappingCached) {
			    std::lock_guard<std::mutex> guard(
				collect.reuseMutex);
			    std::vector<compact_coverage_asset *> &group =
				collect.reuseGroups[broadKey];
			    asset.deferSourceMeshContract =
				asset.sourceMeshRequest.meshAssetName ==
				    asset.dp->d_namep && !group.empty();
			    group.push_back(&asset);
			}
		    }
		}
		bobol_lod_working_set_release(coverageWorkingSet);
	    });

	    if (asset.terminalEmpty) {
		terminalEmptyOccurrences.fetch_add(1);
		continue;
	    }

	    if (asset.coverageReady && !stream->isCancelled()) {
		/*
		 * Publish a synthetic, unselectable draw-target extent before
		 * this asset's leaf box.  Autoview intentionally waits for the
		 * final exact snapshot below, so these monotonic intermediate
		 * extents improve cold visual feedback without moving the camera.
		 */
		SbBox3f overviewSnapshot;
		overviewSnapshot.makeEmpty();
		{
		    std::lock_guard<std::mutex> guard(aggregateMutex);
		    const SbBox3f occurrenceBounds =
			database_source_transform_bounds(
			    asset.coverageBounds,
			    item.occurrence.localTransform);
		    if (!occurrenceBounds.isEmpty())
			aggregateBounds.extendBy(occurrenceBounds);
		    const int64_t now = bu_gettime();
		    if (!aggregateBounds.isEmpty() &&
			(lastOverviewPublication == 0 ||
			 now - lastOverviewPublication >= 100000)) {
			overviewSnapshot = aggregateBounds;
			lastOverviewPublication = now;
		    }
		}
		if (!haveCertifiedWarmBounds && !overviewSnapshot.isEmpty()) {
		    BObolCompactOccurrence overview =
			bobol_database_source_coverage_overview(source, treeName,
			    overviewSnapshot, revision);
		if (overview.geometry)
		    stream->pushPriority(overview);
		}

		const bool includeSourceRequest = asset.sourceMeshReady &&
		    !asset.deferSourceMeshContract;
		BObolCompactOccurrence occurrence =
		    compact_coverage_leaf_occurrence(asset, item.occurrence,
			includeSourceRequest);
		stream->recordManifestOccurrence(occurrence);
		publications.push_back(std::move(occurrence));
		publishedBoxes.fetch_add(1);
		if (asset.deferSourceMeshContract) {
		    std::lock_guard<std::mutex> guard(collect.reuseMutex);
		    asset.deferredOccurrences.push_back(item.occurrence);
		}
	    }

	    /* Only representations which cannot publish a complete lazy contract
	     * from coverage metadata need the deferred realization phase. */
	    if (item.asset && item.asset->dp &&
		item.asset->dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT) {
		std::lock_guard<std::mutex> guard(collect.workMutex);
		collect.detailWork.push_back(std::move(item));
	    }
	    }
	    if (!publications.empty())
		stream->pushBatch(std::move(publications));
	}
    };
    const auto detailWorker = [&]() {
	bu_nice_set(5);
	for (;;) {
	    compact_coverage_work_item item;
	    {
		std::lock_guard<std::mutex> guard(collect.workMutex);
		if (collect.detailWork.empty())
		    break;
		item = std::move(collect.detailWork.front());
		collect.detailWork.pop_front();
	    }
	    if (!item.asset || stream->isCancelled()) {
		stream->notePreparationWorkCompleted();
		continue;
	    }
	    compact_coverage_asset &asset = *item.asset;
	    /* Authored BoTs already published a complete lazy source contract in
	     * the coverage phase.  Their arrays and PoP hierarchy are view-demand
	     * work and must not be imported merely because the leaf exists.  BREP
	     * and analytic primitives still need terminal representation here. */
	    if (asset.dp &&
		asset.dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) {
		stream->notePreparationWorkCompleted();
		continue;
	    }
	    std::call_once(asset.realizeOnce, [&]() {
		if (!bobol_lod_working_set_acquire(
		    asset.estimatedWorkingSetBytes))
		    return;
		struct rt_db_internal intern;
		RT_DB_INTERNAL_INIT(&intern);
		bool ownsInternal = false;
		if (!stream->isCancelled() &&
		    rt_db_get_internal(&intern, asset.dp,
			source->getDatabase(), NULL) >= 0) {
		    ownsInternal = true;
		    const struct rt_bot_internal *bot =
			intern.idb_type == ID_BOT && intern.idb_ptr ?
			static_cast<const struct rt_bot_internal *>(
			    intern.idb_ptr) : NULL;
		    if (bot &&
			bot->num_faces >= source->lodBotThreshold.getValue() &&
			cad_source_mesh_request_from_bot(
			    asset.sourceMeshRequest, bot)) {
			asset.sourceMeshRequest.meshAssetPath =
			    asset.assetPath.c_str();
			asset.sourceMeshRequest.meshAssetName =
			    asset.dp->d_namep ? asset.dp->d_namep : "";
			asset.vertexCount = bot->num_vertices;
			asset.faceCount = bot->num_faces;
			asset.mode = bot->mode;
			asset.orientation = bot->orientation;
			asset.sampleFingerprintValid =
			    compact_stream_lod_sample_fingerprint(
				asset.sampleFingerprint, bot);
			if (asset.coverageReady && asset.coverageGeometry) {
			    asset.geometry = asset.coverageGeometry;
			} else {
			    asset.geometry = bobol_cad_structural_bounds_geometry(
				asset.sourceMeshRequest.bounds,
				asset.proxyGeometryTransform);
			}
			asset.ready = asset.geometry ? true : false;
			asset.sourceMeshReady = asset.ready;
			/* Transfer this already-paid cold import into a bounded,
			 * weakly referenced stream lease.  The first visible LoD
			 * task can build its PoP cache from these arrays instead of
			 * rereading and decoding the same multi-hundred-megabyte
			 * BoT.  If the stream window is disabled or immediately
			 * evicts it, the owner below frees it normally. */
			if (asset.ready && !stream->isCancelled()) {
			    struct rt_db_internal *owned =
				new (std::nothrow) struct rt_db_internal;
			    if (owned) {
				*owned = intern;
				RT_DB_INTERNAL_INIT(&intern);
				ownsInternal = false;
				std::shared_ptr<void> owner(owned,
				    [](void *pointer) {
					struct rt_db_internal *internal =
					    static_cast<struct rt_db_internal *>(
						pointer);
					if (internal) {
					    rt_db_free_internal(internal);
					    delete internal;
					}
				    });
				std::shared_ptr<BObolStagedSourceMesh> staged =
				    std::make_shared<BObolStagedSourceMesh>();
				staged->owner = owner;
				staged->bot =
				    static_cast<const struct rt_bot_internal *>(
					owned->idb_ptr);
				staged->assetName =
				    asset.dp->d_namep ? asset.dp->d_namep : "";
				staged->sourceRevision = revision;
				staged->byteCount =
				    compact_coverage_bot_source_bytes(staged->bot);
				if (stream->retainStagedSource(staged))
				    asset.sourceMeshRequest.stagedSource = staged;
			    }
			}
	    } else if (intern.idb_type == ID_BREP && intern.idb_ptr &&
		source->drawMode.getValue() ==
		    SoBRLDatabaseSource::WIREFRAME) {
		/* BREP wire is a first-class progressive representation.  Its
		 * ranges all reference one immutable line buffer and the renderer
		 * selects the view cut directly.  Building the shaded face-set PoP
		 * as well made one BREP appear twice: correct native curves plus a
		 * coarse triangulation overlay which looked like a lingering box. */
		Obol::PartGeometryBuilder wireGeometry;
		if (cad_progressive_wire_part_geometry_from_provider(
			&intern, source, wireGeometry) > 0 &&
		    wireGeometry.wire) {
		    asset.vertexCount = wireGeometry.wire->segmentCount() * 2;
		    asset.faceCount = 0;
		    asset.proxyGeometryTransform = SbMatrix::identity();
		    asset.geometry = bobol_cad_build_geometry(
			std::move(wireGeometry), "BREP wire asset");
		    asset.ready = asset.geometry ? true : false;
		}
	    } else if (intern.idb_type == ID_BREP && intern.idb_ptr &&
		intern.idb_meth && intern.idb_meth->ft_indexed_face_set) {
		const struct bg_tess_tol ttol = source_tess_tol(source);
		const struct bn_tol tol = BN_TOL_INIT_TOL;
		std::shared_ptr<BObolStagedSourceMesh> staged =
		    cad_staged_mesh_from_primitive_face_set(
			source->getDatabase(), asset.dp, &intern, &ttol, &tol,
			revision, asset.sourceMeshRequest);
		if (staged) {
		    asset.sourceMeshRequest.meshAssetPath =
			asset.assetPath.c_str();
		    asset.sourceMeshRequest.meshAssetName =
			asset.dp->d_namep ? asset.dp->d_namep : "";
		    asset.vertexCount = staged->pointCount;
		    asset.faceCount = staged->faceCount;
		    /* BREP tessellation is already detached and fully owned on this
		     * bounded worker.  Populate the representation-aware PoP cache
		     * here so eviction of the short staged lease cannot strand a
		     * standing box with no database-side BoT fallback. */
		    struct BObolMeshLodCacheStatus cacheStatus =
			BOBOL_MESH_LOD_CACHE_STATUS_INIT;
		    const int stored = bobol_mesh_lod_cache_store_mesh_variant(
			source->getDatabase(), staged->assetName.getString(),
			staged->points, staged->pointCount, staged->normals,
			staged->faces, staged->faceCount, staged->contentKey,
			staged->shadedCullBackfaces, &cacheStatus);
		    if (stored == BRLCAD_OK) {
			asset.sourceMeshRequest.meshAssetContentHash =
			    cacheStatus.cache_key;
			asset.geometry = asset.coverageReady &&
			    asset.coverageGeometry ? asset.coverageGeometry :
			    bobol_cad_structural_bounds_geometry(
				asset.sourceMeshRequest.bounds,
				asset.proxyGeometryTransform);
			asset.ready = asset.geometry ? true : false;
			asset.sourceMeshReady = asset.ready;
		    }
		    if (!asset.ready) {
			/* Tessellation is presentation data; PoP derivation and cache
			 * persistence are optional accelerators.  Preserve a valid
			 * staged mesh if either operation is unavailable or fails.
			 * Discarding it made the compact stream non-authoritative and
			 * forced the fallback hierarchy walk to repeat the same
			 * expensive BREP tessellation before anything but boxes could
			 * be shown. */
			Obol::PartGeometryBuilder terminalGeometry;
			if (cad_mesh_part_geometry_from_staged_source(
				*staged, terminalGeometry)) {
			    asset.geometry = bobol_cad_build_geometry(
				std::move(terminalGeometry),
				"terminal BREP tessellation");
			    asset.ready = asset.geometry ? true : false;
			    asset.sourceMeshReady = false;
			    asset.lodEligible = false;
			}
		    }
		    if (asset.ready && !stream->isCancelled() &&
			stream->retainStagedSource(staged))
			asset.sourceMeshRequest.stagedSource = staged;
		}
	    } else if (intern.idb_ptr) {
		/* Ordinary analytic geometry is terminal at the compact-source
		 * layer.  Its conversion may reach legacy NMG/plot code, so keep
		 * that conversion serialized while the detached workers continue
		 * importing and publishing distinct occurrences concurrently. */
		Obol::PartGeometryBuilder geometry;
		const int internalType = intern.idb_type;
		const int drawMode = source_record_draw_mode(source);
		const bool wireGeometry =
		    primitive_uses_wire_in_mesh_mode(internalType) ||
		    ((drawMode == BOBOL_LOD_DRAW_WIRE ||
		      (drawMode == BOBOL_LOD_DRAW_SHADED_BOTS &&
		       (!intern.idb_meth ||
			!intern.idb_meth->ft_indexed_face_set))) &&
		     internalType != ID_BOT);
		int generated = 0;
		bool viewDependent = false;
		{
		    std::lock_guard<std::mutex> guard(
			collect.primitiveGeometryMutex);
		    if (internalType == ID_PNTS) {
			generated = cad_points_part_geometry_from_pnts(
			    static_cast<const struct rt_pnts_internal *>(
				intern.idb_ptr), geometry);
		    } else if (wireGeometry) {
			generated =
			    cad_wire_part_geometry_from_lod_realization_internal(
				&intern, source, asset.coverageBounds, geometry,
				&viewDependent);
			if (!generated)
			    generated = cad_wire_part_geometry_from_plot_internal(
				&intern, source, geometry);
		    } else {
			generated = cad_mesh_part_geometry_from_internal(
			    &intern, source, geometry);
		    }
		}
		if (generated) {
		    if (drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE)
			(void)cad_mesh_append_hidden_line_edges(geometry);
		    asset.realizedSourceType = primitive_type_label(&intern);
		    asset.realizedGeometryKind = geometry.shaded ? "surface" :
			(geometry.points && !geometry.wire ? "point" : "wire");
		    asset.viewDependentCsgGeometry = viewDependent;
		    asset.vertexCount = geometry.shaded ?
			geometry.shaded->positions.size() :
			(geometry.points ? geometry.points->positions.size() :
			 (geometry.wire ? geometry.wire->segmentPoints.size() : 0));
		    asset.faceCount = geometry.shaded ?
			geometry.shaded->indices.size() / 3 : 0;
		    asset.proxyGeometryTransform = SbMatrix::identity();
		    asset.geometry = bobol_cad_build_geometry(
			std::move(geometry), "direct compact asset");
		    asset.ready = asset.geometry ? true : false;
		}
	    }
		}
		if (ownsInternal)
		    rt_db_free_internal(&intern);
		bobol_lod_working_set_release(
		    asset.estimatedWorkingSetBytes);
	    });

	    if (!asset.ready || stream->isCancelled()) {
		stream->notePreparationWorkCompleted();
		continue;
	    }

	    BObolCompactOccurrence occurrence;
	    occurrence.geometry = asset.geometry;
	    occurrence.geometryTransform =
		asset.proxyGeometryTransform;
	    occurrence.localTransform = item.occurrence.localTransform;
	    occurrence.viewDependentCsgGeometry =
		asset.viewDependentCsgGeometry ? TRUE : FALSE;
	    occurrence.lodBacked = asset.lodEligible ? TRUE : FALSE;
	    /*
	     * The coverage occurrence is not merely a visual AABB: it is the
	     * retained leaf that the view-LoD action will refine in place.  Keep
	     * the source contract discovered by this worker with that occurrence.
	     * Omitting it leaves a leaf marked lodBacked but with no request,
	     * which the submit action must skip.  The completed worker stream is
	     * normally adopted as authoritative, so a later serial realization
	     * cannot be relied upon to repair this metadata.
	     */
	    occurrence.sourceMeshRequestValid = asset.sourceMeshReady ?
		TRUE : FALSE;
	    if (occurrence.sourceMeshRequestValid)
		occurrence.sourceMeshRequest = asset.sourceMeshRequest;
	    occurrence.occurrenceIndex = item.occurrence.occurrenceIndex;
	    occurrence.booleanOperation = item.occurrence.booleanOperation;
	    occurrence.summary = item.occurrence.summary;
	    if (occurrence.sourceMeshRequestValid) {
		compact_source_mesh_request_sync(occurrence.sourceMeshRequest,
		    occurrence.summary);
		compact_summary_lod_from_source_mesh_request(occurrence.summary,
		    occurrence.sourceMeshRequest);
	    } else if (asset.geometry && asset.geometry->wire) {
		occurrence.summary.shapeKind =
		    BObolRealizedShapeSummary::SHAPE_VLIST;
		occurrence.summary.geometryKind =
		    asset.realizedGeometryKind.empty() ? "wire" :
		    asset.realizedGeometryKind.c_str();
		occurrence.summary.pointCount = asset.vertexCount;
		occurrence.summary.commandCount = asset.vertexCount;
		occurrence.summary.segmentCount =
		    asset.geometry->wire->segmentCount();
		occurrence.summary.boundsValid = TRUE;
		occurrence.summary.bounds = asset.geometry->wire->bounds;
	    } else if (asset.geometry && asset.geometry->points) {
		occurrence.summary.shapeKind =
		    BObolRealizedShapeSummary::SHAPE_VLIST;
		occurrence.summary.geometryKind = "point";
		occurrence.summary.pointCount = asset.vertexCount;
		occurrence.summary.commandCount = asset.vertexCount;
		occurrence.summary.pointPrimitiveCount = asset.vertexCount;
		occurrence.summary.boundsValid = TRUE;
		occurrence.summary.bounds = asset.geometry->points->bounds;
	    } else if (asset.geometry && asset.geometry->shaded) {
		occurrence.summary.shapeKind =
		    BObolRealizedShapeSummary::SHAPE_MESH;
		occurrence.summary.geometryKind = "surface";
		occurrence.summary.pointCount = asset.vertexCount;
		occurrence.summary.triangleCount = asset.faceCount;
		occurrence.summary.indexCount =
		    asset.geometry->shaded->indices.size();
		occurrence.summary.boundsValid = TRUE;
		occurrence.summary.bounds = asset.geometry->shaded->bounds;
	    }
	    if (!asset.realizedSourceType.empty())
		occurrence.summary.sourceType =
		    asset.realizedSourceType.c_str();
	    stream->recordManifestOccurrence(occurrence);
	    stream->push(std::move(occurrence));
	    if (!asset.coverageReady)
		publishedTerminalOnly.fetch_add(1);
	    publishedContracts.fetch_add(1);
	    stream->notePreparationWorkCompleted();
	}
    };
    std::vector<std::thread> workers;
    workers.reserve(workerCount);
    if (!selectiveWarmReplay) {
	for (size_t i = 0; i < workerCount; i++)
	    workers.push_back(std::thread(coverageWorker));
    }

    /*
     * Enumeration is the producer, not a prerequisite.  The old two-phase
     * walk held every occurrence until the complete hierarchy was known,
     * creating seconds of blank UI on 50k+ cold scenes while all CPUs sat
     * idle.  Workers now import bounds and publish the expanding overview/
     * leaf boxes while db_walk_tree continues discovering later branches.
     */
    int walkResult = 0;
    if (!selectiveWarmReplay) {
	struct db_tree_state initialState;
	db_init_db_tree_state(&initialState, source->getDatabase());
	initialState.ts_stop_at_regions = 0;
	const char *treeNames[1] = {treeName};
	walkResult = db_walk_tree_leaf_instances(
	    source->getDatabase(), 1, treeNames, 1, &initialState, NULL, NULL,
	    compact_coverage_collect_leaf, &collect);
	db_free_db_tree_state(&initialState);
    }
    /* Accumulate the profile after enumeration, but keep it producer-local
     * until workers distinguish terminal-empty BoTs from invalid nonempty
     * records. Publishing the raw hierarchy count as a complete discovery
     * permanently stranded streams containing valid zero-face placeholders. */
    BObolCompactSourceProfile profile;
    if (walkResult >= 0 && !selectiveWarmReplay) {
	profile.valid = TRUE;
	profile.occurrenceCount = collect.occurrenceCount;
	profile.uniqueAssetCount = collect.assets.size();
	profile.reusedOccurrenceCount = profile.occurrenceCount >
	    profile.uniqueAssetCount ? profile.occurrenceCount -
	    profile.uniqueAssetCount : 0;
	for (const std::unique_ptr<compact_coverage_asset> &asset :
	     collect.assets) {
	    const uint64_t encodedBytes =
		database_source_encoded_source_bytes(source->getDatabase(),
		    asset ? asset->dp : NULL);
	    profile.largestAssetBytes = std::max(profile.largestAssetBytes,
		encodedBytes);
	    profile.encodedSourceBytes = encodedBytes >
		UINT64_MAX - profile.encodedSourceBytes ? UINT64_MAX :
		profile.encodedSourceBytes + encodedBytes;
	}
    }
    {
	std::lock_guard<std::mutex> guard(collect.workMutex);
	collect.work.insert(collect.work.end(),
	    std::make_move_iterator(collect.producerWork.begin()),
	    std::make_move_iterator(collect.producerWork.end()));
	collect.producerWork.clear();
	collect.producerDone = true;
    }
    collect.workReady.notify_all();
    for (std::thread &thread : workers)
	thread.join();
    if (walkResult < 0 || stream->isCancelled())
	return -1;
    const size_t emptyOccurrenceCount = terminalEmptyOccurrences.load();
    const size_t drawableOccurrenceCount =
	emptyOccurrenceCount >= collect.occurrenceCount ? 0 :
	collect.occurrenceCount - emptyOccurrenceCount;
    stream->setExpectedCount(drawableOccurrenceCount);
    if (collect.assets.empty())
	return 0;

    /* Coverage has now produced the exact finite terminal-representation
     * queue.  Publish its denominator before reuse proof or detail workers
     * begin so the owner can distinguish a productive long tessellation from
     * an unranked provider stall. */
    stream->setPreparationWorkCount(collect.detailWork.size());

    /* Resolve only the ambiguous broad groups.  The proof scans borrowed
     * serialized arrays and therefore has a tiny, bounded working set even
     * for multi-gigabyte meshes.  Boxes have already reached the GUI, while
     * only one representative request per plausible group was allowed to
     * enter the expensive resident/PoP pipeline. */
    size_t transformedAssets = 0;
    size_t transformedOccurrences = 0;
    for (auto &groupEntry : collect.reuseGroups) {
	std::vector<compact_coverage_asset *> &group = groupEntry.second;
	if (group.size() < 2)
	    continue;
	std::vector<compact_coverage_asset *> representatives;
	compact_coverage_asset *primary = NULL;
	for (compact_coverage_asset *candidate : group) {
	    if (!candidate)
		continue;
	    if (!candidate->deferSourceMeshContract) {
		representatives.push_back(candidate);
		if (!primary)
		    primary = candidate;
	    }
	}
	if (!primary)
	    continue;

	/* The common expanded-instance case has one representative followed by
	 * many exact copies.  Prove those copies in parallel: each worker scans
	 * read-only mapped pages and owns only small decode blocks.  Candidates
	 * which fail the primary proof are resolved serially below, preserving
	 * correctness for a broad-key collision containing several distinct
	 * duplicate families. */
	std::vector<compact_coverage_asset *> candidates;
	for (compact_coverage_asset *candidate : group) {
	    if (candidate && candidate->deferSourceMeshContract)
		candidates.push_back(candidate);
	}
	std::vector<SbMatrix> primaryTransforms(candidates.size(),
	    SbMatrix::identity());
	std::vector<unsigned char> primaryMatches(candidates.size(), 0);
	std::atomic<size_t> nextCandidate(0);
	const size_t proofWorkerCount = std::max<size_t>(1,
	    std::min(workerCount, candidates.size()));
	std::vector<std::thread> proofWorkers;
	proofWorkers.reserve(proofWorkerCount);
	for (size_t i = 0; i < proofWorkerCount; ++i) {
	    proofWorkers.push_back(std::thread([&]() {
		for (;;) {
		    const size_t candidateIndex =
			nextCandidate.fetch_add(1);
		    if (candidateIndex >= candidates.size())
			break;
		    primaryMatches[candidateIndex] =
			compact_coverage_serialized_rigid_match(
			    source->getDatabase(), primary->dp,
			    candidates[candidateIndex]->dp,
			    primaryTransforms[candidateIndex],
			    mappedDatabase.file) ? 1 : 0;
		}
	    }));
	}
	for (std::thread &proofWorker : proofWorkers)
	    proofWorker.join();

	for (size_t candidateIndex = 0;
	     candidateIndex < candidates.size(); ++candidateIndex) {
	    compact_coverage_asset *candidate = candidates[candidateIndex];
	    const BObolSourceMeshRequest objectRequest =
		candidate->sourceMeshRequest;
	    bool matched = primaryMatches[candidateIndex] != 0;
	    compact_coverage_asset *matchedRepresentative =
		matched ? primary : NULL;
	    SbMatrix transform = primaryTransforms[candidateIndex];
	    if (!matched) {
		for (compact_coverage_asset *representative : representatives) {
		    if (!representative || representative == primary ||
			!representative->sourceMeshReady)
			continue;
		    if (!compact_coverage_serialized_rigid_match(
			    source->getDatabase(), representative->dp,
			    candidate->dp, transform,
			    mappedDatabase.file))
			continue;
		    matchedRepresentative = representative;
		    matched = true;
		    break;
		}
	    }
	    if (matched && matchedRepresentative) {
		candidate->sourceMeshRequest.meshAssetPath =
		    matchedRepresentative->sourceMeshRequest.meshAssetPath;
		candidate->sourceMeshRequest.meshAssetName =
		    matchedRepresentative->sourceMeshRequest.meshAssetName;
		candidate->sourceMeshRequest.meshAssetBounds =
		    matchedRepresentative->sourceMeshRequest.meshAssetBounds;
		candidate->sourceMeshRequest.meshAssetContentHash =
		    matchedRepresentative->sourceMeshRequest.
			meshAssetContentHash;
		candidate->sourceMeshRequest.meshAssetTransform = transform;
		compact_stream_lod_asset_record_store(source->getDatabase(),
		    candidate->dp, objectRequest,
		    matchedRepresentative->sourceMeshRequest, transform);
		transformedAssets++;
	    }
	    if (!matched)
		representatives.push_back(candidate);

	    /* Whether reused or independently retained, release the source
	     * contract only after the ambiguity has been resolved.  Merging this
	     * richer record preserves the already-visible box, placement,
	     * selection state and semantic identity. */
	    candidate->deferSourceMeshContract = false;
	    for (const compact_coverage_occurrence &leaf :
		 candidate->deferredOccurrences) {
		BObolCompactOccurrence upgrade =
		    compact_coverage_leaf_occurrence(*candidate, leaf, true);
		stream->recordManifestOccurrence(upgrade);
		stream->push(std::move(upgrade));
		publishedContracts.fetch_add(1);
		transformedOccurrences++;
	    }
	    candidate->deferredOccurrences.clear();
	}
    }
    if (getenv("BOBOL_DRAW_TIMING") && transformedOccurrences)
	bu_log("[obol-timing] serialized transformed reuse: %zu assets, "
	       "%zu contracts\n", transformedAssets,
	       transformedOccurrences);

    /*
     * Leaf AABBs may already fill the producer queue faster than the GUI can
     * adopt them.  Once the parallel bound pass knows every occurrence,
     * publish one unselectable whole-target extent through the priority lane.
     * It is intentionally absent from the final authoritative compact index;
     * adoption therefore removes it atomically after all leaf boxes/meshes
     * are present.
     */
    if (!aggregateBounds.isEmpty()) {
	/* The aggregate overview remains useful even when it is conservative:
	 * it is the earliest complete visual scope of the draw target.  A warm
	 * exact extent is stronger evidence and must never be replaced by this
	 * provisional union. */
	if (!haveCertifiedWarmBounds) {
	    stream->setCoverageBounds(aggregateBounds);
	    BObolCompactOccurrence overview =
		bobol_database_source_coverage_overview(source, treeName,
		    aggregateBounds, revision);
	    if (overview.geometry)
		stream->pushPriority(overview);
	}

	SbBox3f terminalBounds = haveCertifiedWarmBounds ?
	    certifiedWarmBounds : aggregateBounds;
	bool terminalBoundsExact = haveCertifiedWarmBounds ||
	    collect.aggregateBoundsExact;
	if (!terminalBoundsExact) {
	    SbBox3f semanticBounds;
	    if (source_bounds_from_database_path(source, treeName,
		    semanticBounds)) {
		terminalBounds = semanticBounds;
		terminalBoundsExact = true;
	    }
	}
	if (terminalBoundsExact) {
	    stream->setCoverageBounds(terminalBounds);
	    if (!haveCertifiedWarmBounds && terminalBounds != aggregateBounds) {
		BObolCompactOccurrence exactOverview =
		    bobol_database_source_coverage_overview(source, treeName,
			terminalBounds, revision);
		if (exactOverview.geometry)
		    stream->pushPriority(exactOverview);
	}
	    /* Completion means semantic exactness, not merely that every leaf
	     * contributed to a conservative aggregate.  The priority lane is
	     * drained first, so a consumer observing this flag will publish the
	     * terminal overview before applying the one-shot autoview. */
	    stream->setCoverageBoundsComplete(true);
	}
    }
    const int64_t coverageCompleted = bu_gettime();
    if (getenv("BOBOL_DRAW_TIMING"))
	bu_log("[obol-timing] coverage bounds: %.1f ms; %zu occurrences, "
	       "%zu assets, %zu workers, %zu boxes\n",
	       static_cast<double>(coverageCompleted - collectStart) / 1000.0,
	       collect.occurrenceCount, collect.assets.size(), workerCount,
	       publishedBoxes.load());

    /* Bounds-first is a scheduling barrier, not a loss of parallelism.  Once
     * complete visual coverage and an exact target extent are available,
     * reuse the same bounded worker count to import the source arrays needed
     * by PoP and enrich each standing leaf in place. */
    workers.clear();
    for (size_t i = 0; i < workerCount; i++)
	workers.push_back(std::thread(detailWorker));
    for (std::thread &thread : workers)
	thread.join();
    if (stream->isCancelled())
	return -1;
    if (profile.valid) {
	profile.occurrenceCount = drawableOccurrenceCount;
	profile.uniqueAssetCount = 0;
	profile.encodedSourceBytes = 0;
	profile.largestAssetBytes = 0;
	for (const std::unique_ptr<compact_coverage_asset> &asset :
	     collect.assets) {
	    if (!asset || asset->terminalEmpty)
		continue;
	    profile.uniqueAssetCount++;
	    const uint64_t encodedBytes =
		database_source_encoded_source_bytes(source->getDatabase(),
		    asset->dp);
	    profile.largestAssetBytes = std::max(
		profile.largestAssetBytes, encodedBytes);
	    profile.encodedSourceBytes = encodedBytes >
		UINT64_MAX - profile.encodedSourceBytes ? UINT64_MAX :
		profile.encodedSourceBytes + encodedBytes;
	}
	profile.reusedOccurrenceCount = profile.occurrenceCount >
	    profile.uniqueAssetCount ? profile.occurrenceCount -
	    profile.uniqueAssetCount : 0;
	profile.valid = profile.occurrenceCount > 0 &&
	    profile.uniqueAssetCount > 0 &&
	    profile.uniqueAssetCount <= profile.occurrenceCount &&
	    profile.encodedSourceBytes > 0 &&
	    profile.largestAssetBytes > 0;
	if (profile.isValid(drawableOccurrenceCount)) {
	    stream->setSourceProfile(profile);
	}
    }
    const int64_t detailCompleted = bu_gettime();

    /* A unique broad signature cannot be a rigid transformed copy of another
     * asset in this root.  Seed those self-assets into the realization cache,
     * allowing the authoritative serial identity walk below to avoid a second
     * import.  Ambiguous groups retain only their already-published boxes and
     * proceed through exact PCA verification before sharing geometry. */
    std::unordered_map<std::string, size_t> broadCounts;
    for (const std::unique_ptr<compact_coverage_asset> &asset :
	 collect.assets) {
	const std::string key = compact_coverage_broad_key(*asset);
	if (!key.empty())
	    broadCounts[key]++;
    }
    size_t seeded = 0;
    for (const std::unique_ptr<compact_coverage_asset> &assetPtr :
	 collect.assets) {
	compact_coverage_asset &asset = *assetPtr;
	const std::string broadKey = compact_coverage_broad_key(asset);
	const bool brepAsset = asset.dp &&
	    asset.dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP;
	const bool nativeBrepWire = brepAsset && asset.geometry &&
	    asset.geometry->wire && !asset.sourceMeshReady;
	if ((!brepAsset &&
	     (broadKey.empty() || broadCounts[broadKey] != 1)) ||
	    !asset.geometry || !asset.ready ||
	    (!asset.sourceMeshReady && !nativeBrepWire))
	    continue;
	/* A self asset needs no transformed-reuse mapping record.  The completed
	 * leaf manifest persists its canonical source request directly; writing
	 * one LMDB record for every ordinary unique mesh made cold completion
	 * transaction-bound and inflated the draw cache by tens of megabytes.
	 *
	 * The cached geometry is the shared unit AABB used by cold coverage, not
	 * an object-coordinate box.  Preserve its geometry-to-asset transform.
	 * Storing identity here made the following authoritative walk install a
	 * 0..1 proxy while retaining the real mesh request bounds.  Autoview then
	 * framed the unit cube and view culling correctly rejected the actual
	 * (often million-unit) mesh forever.
	 */
	const SbBox3f &cachedBounds = asset.sourceMeshReady ?
	    asset.sourceMeshRequest.bounds : asset.geometry->wire->bounds;
	cache->storeMeshCadGeometryReference(asset.cacheKey, asset.geometry,
	    asset.proxyGeometryTransform, brepAsset ? "brep" : "bot",
	    nativeBrepWire ? "wire" : "aabb", &cachedBounds, true,
	    asset.sourceMeshReady ? &asset.sourceMeshRequest : NULL);
	seeded++;
    }
    if (getenv("BOBOL_DRAW_TIMING"))
	bu_log("[obol-timing] coverage prepass: %.1f ms "
	       "(bounds %.1f, detail %.1f); %zu occurrences, "
	       "%zu assets, %zu workers, %zu boxes, %zu contracts, "
	       "%zu unique assets seeded, "
	       "global peak=%zu bytes/%zu tasks\n",
	       static_cast<double>(bu_gettime() - collectStart) / 1000.0,
	       static_cast<double>(coverageCompleted - collectStart) / 1000.0,
	       static_cast<double>(detailCompleted - coverageCompleted) / 1000.0,
	       collect.occurrenceCount, collect.assets.size(), workerCount,
	       publishedBoxes.load(), publishedContracts.load(),
	       seeded,
	       bobol_lod_working_set_global_peak_bytes(),
	       bobol_lod_working_set_global_peak_tasks());

    /* Every occurrence now has either a lazy source-mesh contract or terminal
     * compact geometry.  When all assets reached one of those states, this
     * stream is the authoritative realization and a second serial hierarchy
     * walk would only repeat imports, tessellation, identity construction and
     * publication already completed above. */
    const size_t publishedOccurrenceCoverage =
	publishedBoxes.load() + publishedTerminalOnly.load();
    bool authoritativeStream = drawableOccurrenceCount > 0 &&
	publishedOccurrenceCoverage == drawableOccurrenceCount;
    size_t uncertifiedCoverage = 0;
    size_t uncertifiedRepresentation = 0;
    for (const std::unique_ptr<compact_coverage_asset> &asset :
	 collect.assets) {
	if (!asset || !asset->dp) {
	    uncertifiedCoverage++;
	    authoritativeStream = false;
	    continue;
	}
	if (asset->terminalEmpty)
	    continue;
	const bool lazyMesh = asset->lodEligible && asset->sourceMeshReady;
	const bool terminalGeometry = asset->ready && asset->geometry &&
	    !asset->geometry->structuralProxy;
	if (!asset->coverageReady && !terminalGeometry) {
	    uncertifiedCoverage++;
	    if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
		bu_log("[obol-timing] uncertified coverage asset: %s "
		       "minor=%d ready=%d\n",
		       asset->dp->d_namep ? asset->dp->d_namep : "?",
		       asset->dp->d_minor_type,
		       asset->coverageReady ? 1 : 0);
	    authoritativeStream = false;
	    continue;
	}
	if (!lazyMesh && !terminalGeometry) {
	    uncertifiedRepresentation++;
	    if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
		bu_log("[obol-timing] uncertified representation asset: %s "
		       "minor=%d lod=%d source=%d ready=%d geometry=%d "
		       "structural=%d\n",
		       asset->dp->d_namep ? asset->dp->d_namep : "?",
		       asset->dp->d_minor_type,
		       asset->lodEligible ? 1 : 0,
		       asset->sourceMeshReady ? 1 : 0,
		       asset->ready ? 1 : 0,
		       asset->geometry ? 1 : 0,
		       asset->geometry && asset->geometry->structuralProxy ? 1 : 0);
	    authoritativeStream = false;
	}
    }
    if (getenv("BOBOL_DRAW_TIMING"))
	bu_log("[obol-timing] coverage certification: authoritative=%d "
	       "occurrences=%zu drawable=%zu empty=%zu boxes=%zu "
	       "terminal_only=%zu assets=%zu "
	       "missing_bounds=%zu "
	       "missing_terminal=%zu\n", authoritativeStream ? 1 : 0,
	       collect.occurrenceCount, drawableOccurrenceCount,
	       emptyOccurrenceCount, publishedBoxes.load(),
	       publishedTerminalOnly.load(), collect.assets.size(),
	       uncertifiedCoverage,
	       uncertifiedRepresentation);
    if (authoritativeStreamOut)
	*authoritativeStreamOut = authoritativeStream;
    (void)stream->sealManifest(drawableOccurrenceCount);
    return publishedBoxes.load() > static_cast<size_t>(INT_MAX) ? INT_MAX :
	static_cast<int>(publishedBoxes.load());
}

} /* anonymous namespace */

BObolCompactOccurrence
bobol_database_source_coverage_overview(const SoBRLDatabaseSource *source,
	const char *treeName, const SbBox3f &bounds, uint32_t revision)
{
    BObolCompactOccurrence overview;
    if (!source || !treeName || !treeName[0] || bounds.isEmpty())
	return overview;

    overview.geometry = compact_coverage_overview_geometry(bounds);
    overview.geometryTransform = SbMatrix::identity();
    if (!overview.geometry)
	return overview;
    overview.summary = compact_occurrence_summary(source,
	treeName, treeName, "proxy", "overview-aabb", revision,
	BObolRealizedShapeSummary::SHAPE_VLIST);
    overview.summary.recordRole = "lod-overview";
    overview.summary.visible = TRUE;
    overview.summary.selectable = FALSE;
    overview.summary.lodAvailable = TRUE;
    overview.summary.lodActiveCut = BOBOL_LOD_QUALITY_PROXY;
    overview.summary.lodBoundsMin = bounds.getMin();
    overview.summary.lodBoundsMax = bounds.getMax();
    overview.summary.pointCount = 24;
    overview.summary.commandCount = 24;
    overview.summary.segmentCount = 12;
    overview.summary.boundsValid = TRUE;
    overview.summary.bounds = bounds;
    return overview;
}

bool
bobol_database_source_uses_serialized_bot_coverage(
    const SoBRLDatabaseSource *source, struct db_i *database)
{
    if (!source || !database || db_version(database) != 5 ||
	!source->usesMeshRealization() || source->lodBotThreshold.getValue() == 0 ||
	source_uses_evaluated_path_realization(source) ||
	source_has_auxiliary_children(source))
	return false;

    const std::string name = database_source_leaf_component(source->path.getValue());
    struct directory *dp = name.empty() ? RT_DIR_NULL :
	db_lookup(database, name.c_str(), LOOKUP_QUIET);
    return dp && !(dp->d_flags & RT_DIR_COMB) &&
	dp->d_major_type == DB5_MAJORTYPE_BRLCAD &&
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT;
}

int
bobol_database_source_construct_mesh_compact_with_cache(
    SoBRLDatabaseSource *source,
    BObolDatabaseSourceRealizationCache *cache,
    BObolCompactOccurrenceStream *stream)
{
    BObolPerformanceTimer timer(BOBOL_PERF_MESH_REALIZE_US);
    if (timer.active())
	bobol_performance_counter_add(BOBOL_PERF_MESH_REALIZE_CALLS, 1);

    BObolDatabaseSourceRealizationCache localCache;
    if (!cache)
	cache = &localCache;

    if (!source)
	return -1;

    if (source_uses_evaluated_path_realization(source))
	return 0;

    source->d->compactHandleSourceId = source_stable_compact_handle_id(source);

    source->realizationDiagnostic = "";
    if (!source->d->dbip) {
	source->realizationDiagnostic = "database source has no database";
	return -1;
    }

    std::string treeNameStorage =
	database_lookup_path_from_source_path(source->path.getValue());
    const char *treeName = treeNameStorage.c_str();
    if (!treeName[0]) {
	source->realizationDiagnostic = "database source path is empty";
	return -1;
    }

    std::shared_ptr<std::mutex> lodRealizationMutex;
    std::unique_lock<std::mutex> lodRealizationLock;
    if (source->lodBotThreshold.getValue() > 0) {
	lodRealizationMutex = compact_stream_lod_realization_mutex(
	    source->d->dbip, treeName);
	const int64_t waitStart = bu_gettime();
	lodRealizationLock = std::unique_lock<std::mutex>(
	    *lodRealizationMutex);
	const int64_t waitMicroseconds = bu_gettime() - waitStart;
	if (getenv("BOBOL_DRAW_TIMING") && waitMicroseconds >= 1000)
	    bu_log("[obol-timing] stream LoD single-flight: root=%s "
		   "stream=%d waited %.1f ms\n", treeName, stream ? 1 : 0,
		   static_cast<double>(waitMicroseconds) / 1000.0);
    }

    (void)treeName;
    const uint32_t revision = source->sourceRevision.getValue();
    const bool serializedBotCoverage = stream &&
	bobol_database_source_uses_serialized_bot_coverage(source, source->d->dbip);
    int directRealized = 0;
    if (!serializedBotCoverage) {
	BObolPerformanceTimer directTimer(BOBOL_PERF_DIRECT_LEAF_US);
	if (directTimer.active())
	    bobol_performance_counter_add(BOBOL_PERF_DIRECT_LEAF_CALLS, 1);
	directRealized = realize_direct_leaf_mesh_compact(source, cache,
	    revision, stream);
	if (directTimer.active()) {
	    if (directRealized > 0) {
		bobol_performance_counter_add(
		    BOBOL_PERF_DIRECT_LEAF_REALIZED, 1);
	    } else if (directRealized < 0) {
		bobol_performance_counter_add(
		    BOBOL_PERF_DIRECT_LEAF_FAILED, 1);
	    } else {
		bobol_performance_counter_add(
		    BOBOL_PERF_DIRECT_LEAF_FALLBACK, 1);
	    }
	}
    }

    if (directRealized > 0) {
	SbBox3f semanticBounds;
	if (source_bounds_from_database_path(source, treeName, semanticBounds))
	    (void)source->setSourceBoundsState(TRUE, semanticBounds.getMin(),
		semanticBounds.getMax(), TRUE);
	mark_source_realized_current(source);
	return 1;
    }
    if (directRealized < 0) {
	remove_non_auxiliary_children(source);
	source->discardCompactInstanceHistory();
	source->realizationIdentity = "";
	return -1;
    }

    /* Prefill is a batch optimization for synchronous static realization.
     * It is the wrong latency and memory tradeoff for progressive BoTs: the
     * occurrence walk already performs bounded online PCA while structural
     * boxes remain visible.  Running both paths can import two candidate/
     * representative pairs concurrently. */
    if (!stream && source->lodBotThreshold.getValue() == 0 &&
	compact_mesh_prefill_cache(source, cache, treeName, NULL) < 0) {
	remove_non_auxiliary_children(source);
	source->discardCompactInstanceHistory();
	source->realizationIdentity = "";
	return -1;
    }

    bool authoritativeStream = false;
    int streamedCoverage = 0;
    if (stream && source->lodBotThreshold.getValue() > 0 &&
	!stream->hasWarmCoverageComplete()) {
	streamedCoverage = compact_stream_publish_parallel_coverage(source,
	    cache, treeName, stream, revision, &authoritativeStream);
    }
    if (streamedCoverage < 0) {
	remove_non_auxiliary_children(source);
	source->discardCompactInstanceHistory();
	source->realizationIdentity = "";
	return -1;
    }

    if (authoritativeStream) {
	/* The detached worker has no reason to retain a second copy of the
	 * streamed compact registry.  Publish only its terminal source contract;
	 * adoptDetachedCompactRealization(..., TRUE) will certify the drained live
	 * registry and retire the temporary whole-target overview in place. */
	SbBox3f bounds;
	if (stream->getCoverageBounds(bounds) && !bounds.isEmpty()) {
	    (void)source->setSourceBoundsState(TRUE, bounds.getMin(),
		bounds.getMax(), TRUE);
	}
	if (stream->hasWarmCensusComplete())
	    stream->setWarmCoverageComplete(true);
	mark_source_realized_current(source);
	bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_SOURCES, 1);
	bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_INSTANCES,
	    static_cast<uint64_t>(streamedCoverage));
	return streamedCoverage;
    }

    if (serializedBotCoverage) {
	/* Source admission reserved directory/coverage work, not an unrestricted
	 * rt_db_get_internal import.  A malformed or incomplete serialized census
	 * must fail at this boundary instead of entering the serial fallback. */
	source->realizationDiagnostic =
	    "serialized BoT coverage did not produce a complete representation";
	return -1;
    }

    struct db_tree_state init_state;
    db_init_db_tree_state(&init_state, source->d->dbip);
    init_state.ts_stop_at_regions = 0;

    BObolMaterialColorSweep materialSweep(source->d->dbip);
    realize_walk_data data;
    data.source = source;
    data.cache = cache;
    data.revision = revision;
    data.compact_index = std::make_unique<BObolCompactInstanceIndex>();
    data.stream_sink = stream;
    data.material_sweep = &materialSweep;

    const char *av[1] = {treeName};
    const int ret = db_walk_tree_leaf_instances(source->d->dbip, 1, av, 1,
	&init_state, NULL, NULL, realize_mesh_leaf, &data);
    db_free_db_tree_state(&init_state);
    if (getenv("BOBOL_DRAW_TIMING") &&
	source->lodBotThreshold.getValue() > 0)
	bu_log("[obol-timing] stream LoD reuse: cache=%zu deferred=%zu "
	       "pca=%zu reused=%zu representatives=%zu representative-imports=%zu\n",
	       data.stream_lod_asset_hits, data.stream_lod_pca_deferred,
	       data.stream_lod_pca_evaluated, data.stream_lod_pca_reused,
	       data.stream_lod_reuse.size(),
	       data.stream_lod_representative_imports);

    if (ret < 0 || data.realized_shapes <= 0 || data.failed_shapes > 0 ||
	data.compact_unsupported || data.compact_index->entries.empty()) {
	remove_non_auxiliary_children(source);
	source->discardCompactInstanceHistory();
	source->realizationIdentity = "";
	if (data.diagnostic.getLength() > 0)
	    source->realizationDiagnostic = data.diagnostic;
	else {
	    SbString diagnostic;
	    diagnostic.sprintf(
		"%s: compact mesh realization produced no usable occurrences",
		treeName);
	    source->realizationDiagnostic = diagnostic;
	}
	return ret < 0 || data.failed_shapes > 0 ? -1 : 0;
    }

    source->installCompactInstanceIndex(std::move(data.compact_index), TRUE,
	source->d->compactHandleSourceId);
    source->markCompiledAssemblyDirty();

    if (stream)
	(void)close_compact_stream_contract(source, stream);

    const SbBool preserveExact = source->hasExactSourceBounds();

    SbBox3f semanticBounds;
    if (!preserveExact && source_bounds_from_database_path(source, treeName,
	    semanticBounds))
	(void)source->setSourceBoundsState(TRUE, semanticBounds.getMin(),
	    semanticBounds.getMax(), TRUE);
    else if (data.compact_bounds_valid && !data.compact_bounds.isEmpty() &&
	!preserveExact)
	set_source_bounds_from_presentation(source, data.compact_bounds);
    else if (!preserveExact &&
	(!data.compact_bounds_valid || data.compact_bounds.isEmpty()))
	source->clearSourceBounds();
    mark_source_realized_current(source);
    bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_SOURCES, 1);
    bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_INSTANCES,
	static_cast<uint64_t>(source->d->compactIndex->entries.size()));
    return 1;
}

SbBool
bobol_database_source_construct_mesh_with_cache(
    SoBRLDatabaseSource *source,
    BObolDatabaseSourceRealizationCache *cache)
{
    BObolPerformanceTimer timer(BOBOL_PERF_MESH_REALIZE_US);
    if (timer.active())
	bobol_performance_counter_add(BOBOL_PERF_MESH_REALIZE_CALLS, 1);

    BObolDatabaseSourceRealizationCache localCache;
    if (!cache)
	cache = &localCache;

    if (!source)
	return FALSE;

    source->realizationDiagnostic = "";
    if (!source->d->dbip) {
	source->realizationDiagnostic = "database source has no database";
	return FALSE;
    }

    std::string treeNameStorage =
	database_lookup_path_from_source_path(source->path.getValue());
    const char *treeName = treeNameStorage.c_str();
    if (!treeName[0]) {
	source->realizationDiagnostic = "database source path is empty";
	return FALSE;
    }

    remove_non_auxiliary_children(source);
    source->discardCompactInstanceHistory();
    (void)remove_source_placement_transform(source);

    const uint32_t revision = source->sourceRevision.getValue();
    if (source_uses_evaluated_points_realization(source)) {
	if (realize_evaluated_points_source(source, revision) > 0) {
	    mark_source_realized_current(source);
	    SbBox3f sourceBounds;
	    if (!source->getSourceBounds(sourceBounds))
		update_source_bounds_from_realized_geometry(source, FALSE);
	    source->syncRealizedShapeOwnerState();
	    return TRUE;
	}
	remove_non_auxiliary_children(source);
	source->realizationIdentity = "";
	return FALSE;
    }

    int directRealized = 0;
    {
	BObolPerformanceTimer directTimer(BOBOL_PERF_DIRECT_LEAF_US);
	if (directTimer.active())
	    bobol_performance_counter_add(BOBOL_PERF_DIRECT_LEAF_CALLS, 1);
	directRealized = realize_direct_leaf_mesh(source, cache, revision);
	if (directTimer.active()) {
	    if (directRealized > 0) {
		bobol_performance_counter_add(
		    BOBOL_PERF_DIRECT_LEAF_REALIZED, 1);
	    } else if (directRealized < 0) {
		bobol_performance_counter_add(
		    BOBOL_PERF_DIRECT_LEAF_FAILED, 1);
	    } else {
		bobol_performance_counter_add(
		    BOBOL_PERF_DIRECT_LEAF_FALLBACK, 1);
	    }
	}
    }
    if (directRealized > 0) {
	SbBox3f semanticBounds;
	if (source_bounds_from_database_path(source, treeName, semanticBounds))
	    (void)source->setSourceBoundsState(TRUE, semanticBounds.getMin(),
		semanticBounds.getMax(), TRUE);
	mark_source_realized_current(source);
	update_source_bounds_from_realized_geometry(source);
	source->syncRealizedShapeOwnerState();
	return TRUE;
    }
    if (directRealized < 0) {
	remove_non_auxiliary_children(source);
	source->realizationIdentity = "";
	return FALSE;
    }

    struct db_tree_state init_state;
    db_init_db_tree_state(&init_state, source->d->dbip);
    init_state.ts_stop_at_regions = 0;

    struct realize_walk_data data;
    data.source = source;
    data.cache = cache;
    data.revision = source->sourceRevision.getValue();
    data.realized_shapes = 0;
    data.failed_shapes = 0;

    const char *av[1] = { treeName };
    int ret = db_walk_tree_leaf_instances(source->d->dbip, 1, av, 1, &init_state,
					  NULL, NULL, realize_mesh_leaf, &data);
    db_free_db_tree_state(&init_state);

    if (ret < 0 || data.realized_shapes <= 0 || data.failed_shapes > 0) {
	remove_non_auxiliary_children(source);
	source->realizationIdentity = "";
	if (data.diagnostic.getLength() > 0) {
	    source->realizationDiagnostic = data.diagnostic;
	} else if (data.realized_shapes <= 0) {
	    SbString msg;
	    msg.sprintf("%s: no drawable mesh geometry realized", treeName);
	    source->realizationDiagnostic = msg;
	} else {
	    SbString msg;
	    msg.sprintf("%s: mesh realization failed", treeName);
	    source->realizationDiagnostic = msg;
	}
	return FALSE;
    }

    mark_source_realized_current(source);
    update_source_bounds_from_realized_geometry(source);
    source->syncRealizedShapeOwnerState();
    return TRUE;
}
static int
construct_prototype_wireframe(SoBRLDatabaseSource *source)
{
    SbModernUtils::SoNodeRef owner(new SoBRLVListShape);
    auto *shape = static_cast<SoBRLVListShape *>(owner.get());
    /* A revision-dependent square makes an obsolete diagnostic easy to see. */
    constexpr uint32_t extentCycle = 4;
    constexpr float minimumHalfExtent = 1.0f;
    constexpr float halfExtentStep = 0.25f;
    constexpr int pointCount = 5;
    const float halfExtent = minimumHalfExtent + halfExtentStep *
	static_cast<float>(source->sourceRevision.getValue() % extentCycle);
    const SbVec3f points[pointCount] = {
	SbVec3f(-halfExtent, -halfExtent, 0.0f),
	SbVec3f(halfExtent, -halfExtent, 0.0f),
	SbVec3f(halfExtent,  halfExtent, 0.0f),
	SbVec3f(-halfExtent,  halfExtent, 0.0f),
	SbVec3f(-halfExtent, -halfExtent, 0.0f)
    };
    const int32_t commands[pointCount] = {
	SoBRLVListShape::MOVE,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW
    };

    shape->sourcePath = source->path.getValue();
    const char *leafName = lookup_name_from_path(source->path.getValue());
    shape->sourceName = leafName ? leafName : source->path.getValue().getString();
    shape->sourceType = "prototype";
    shape->sourceId = source->sourceRevision.getValue();
    shape->displayName = source->displayName.getValue().getLength() > 0 ?
			 source->displayName.getValue() : shape->sourceName.getValue();
    shape->geometryName = shape->sourceName.getValue();
    shape->sourceIdentity = source_record_identity(source,
			    shape->sourcePath.getValue().getString());
    shape->cacheIdentity = record_identity_with_revision(
			       shape->sourceIdentity.getValue().getString(),
			       shape->sourceId.getValue());
    shape->databaseIntent = FALSE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = TRUE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = BOBOL_LOD_DRAW_DIAGNOSTIC;
    shape->recordRole = "prototype";
    shape->geometryKind = "";
    shape->regionId = 0;
    shape->airCode = 0;
    shape->materialId = 0;
    shape->los = 0;
    shape->materialColorValid = FALSE;
    shape->materialColor = SbColor(1.0f, 1.0f, 1.0f);
    shape->materialShader = "";
    shape->setLineSet(points, commands, pointCount);
    database_source_add_realized_child(source, shape);
    (void)source->setSourceBoundsState(TRUE,
	SbVec3f(-halfExtent, -halfExtent, 0.0f),
	SbVec3f(halfExtent, halfExtent, 0.0f), TRUE);
    mark_source_realized_current(source);
    sync_realized_shape_owner_state_in_node(shape, source);
    return TRUE;
}

static const char *
external_string_or_default(const char *value, const char *fallback)
{
    return value && value[0] ? value : fallback;
}

static const char *
external_source_leaf_name(const SoBRLDatabaseSource *source)
{
    if (!source)
	return "";

    const char *sourcePath = source->path.getValue().getString();
    const char *leaf = lookup_name_from_path(source->path.getValue());
    if (leaf && leaf[0])
	return leaf;
    return sourcePath ? sourcePath : "";
}

template <typename ShapeT>
static void
assign_external_primary_identity(ShapeT *shape,
				 const SoBRLDatabaseSource *source,
				 const char *sourceType,
				 const char *geometryKind)
{
    if (!shape || !source)
	return;

    const char *sourcePath = source->path.getValue().getString();
    const char *sourceName = external_source_leaf_name(source);
    assign_realized_identity(shape, NULL, sourcePath, sourceName, sourceType,
			     source->sourceRevision.getValue(), source);
    shape->geometryKind = geometryKind ? geometryKind : "";
    sync_shape_placement_state(shape, source);
}

static int
external_vlist_command_valid(int32_t command)
{
    return command == SoBRLVListShape::MOVE ||
	   command == SoBRLVListShape::DRAW ||
	   command == SoBRLVListShape::POINT;
}

static SbBool
external_bounds_from_points(const SbVec3f *points,
			    int count,
			    SbVec3f &boundsMin,
			    SbVec3f &boundsMax)
{
    if (!points || count <= 0)
	return FALSE;

    boundsMin = points[0];
    boundsMax = points[0];
    for (int i = 1; i < count; i++) {
	for (int axis = 0; axis < 3; axis++) {
	    if (points[i][axis] < boundsMin[axis])
		boundsMin[axis] = points[i][axis];
	    if (points[i][axis] > boundsMax[axis])
		boundsMax[axis] = points[i][axis];
	}
    }

    return TRUE;
}

template <typename Shape>
static void
prepare_realized_shape_owner(Shape &shape, const SoBRLDatabaseSource &source,
    uint32_t revision)
{
    sync_shape_source_state(&shape, &source);
    shape.ownerSourceRevision = revision;
    shape.ownerRealizationIdentity = source_realization_identity(&source, revision);
    const SourceRealizationValues values{{{revision,
	revision, source.inputsRevision.getValue(),
	source.viewRevision.getValue()}}, SoBRLDatabaseSource::REALIZED, "",
	FALSE, SoBRLDatabaseSource::STALE_NONE};
    PreparedShapeRealization publication(shape, values);
    publication.fields.commit();
}

template <typename Shape>
static void
prepare_realized_shape_owner(Shape &shape, const SoBRLDatabaseSource &source)
{
    prepare_realized_shape_owner(shape, source, source.sourceRevision.getValue());
}

static void
clear_external_shape_geometry(SoBRLVListShape &shape)
{
    shape.setLineSet(NULL, NULL, 0);
    shape.setPrecisePoints(NULL, 0);
}

static void
clear_external_shape_geometry(SoBRLMeshShape &shape)
{
    shape.setIndexedTriangles(NULL, 0, NULL, 0);
}

static std::vector<SoField *>
external_geometry_fields(SoBRLVListShape &shape)
{
    return {&shape.point, &shape.command, &shape.annotationPoint,
	&shape.annotationSegmentTextValid, &shape.annotationSegmentKind,
	&shape.annotationSegmentStart, &shape.annotationSegmentEnd,
	&shape.annotationTextRefPoint, &shape.annotationText,
	&shape.pointColorValid, &shape.pointColor, &shape.pointScaleValid,
	&shape.pointScale, &shape.pointNormalValid, &shape.pointNormal};
}

static std::vector<SoField *>
external_geometry_fields(SoBRLMeshShape &shape)
{
    return {&shape.point, &shape.coordIndex, &shape.normal};
}

enum class OwnedShapeState { Realized, Configured, DisplayName, MaterialPolicy, Placement, Owner };

/* Keep the public patch's field mapping in one place for filtering source
 * no-ops and applying the changed properties to detached source/shape fields. */
template <typename Node, typename Patch, typename Visit>
static void
visit_display_patch(Node &node, Patch &patch, Visit visit)
{
    visit(node.visible, patch.visibleValid, patch.visible);
    visit(node.selected, patch.selectedValid, patch.selected);
    visit(node.highlighted, patch.highlightedValid, patch.highlighted);
    visit(node.lineStyle, patch.lineStyleValid, patch.lineStyle);
    visit(node.lineWidth, patch.lineWidthValid, patch.lineWidth);
    visit(node.transparency, patch.transparencyValid, patch.transparency);
    visit(node.colorOverride, patch.colorOverrideValid, patch.colorOverride);
    visit(node.color, patch.colorValid, patch.color);
    visit(node.selectedColor, patch.selectedColorValid, patch.selectedColor);
    visit(node.highlightedColor, patch.highlightedColorValid, patch.highlightedColor);
    visit(node.ghostedColor, patch.ghostedColorValid, patch.ghostedColor);
}

template <typename T>
static bool publication_value_different(const T &a, const T &b) { return a != b; }
static bool publication_value_different(float a, float b) { return database_source_float_different(a, b); }
static bool publication_value_different(const SbColor &a, const SbColor &b) { return !database_source_color_equal(a, b); }

static bool
filter_display_patch(const SoBRLDatabaseSource &source, BObolDatabaseSourceDisplayPatch &patch)
{
    bool changed = false;
    visit_display_patch(source, patch, [&](const auto &field, SbBool &valid, const auto &value) {
	valid = valid && publication_value_different(field.getValue(), value);
	changed = changed || valid;
    });
    return changed;
}

template <typename Node>
static void
apply_display_patch(Node &node, const BObolDatabaseSourceDisplayPatch &patch)
{
    visit_display_patch(node, patch, [](auto &field, SbBool valid, const auto &value) {
	if (valid) field.setValue(value);
    });
}

struct SourceDisplayChanges {
    BObolDatabaseSourceDisplayPatch patch;
    const BObolDatabaseSourcePresentationPatch *presentation = nullptr;
    bool materialColor = false;
    bool materialRevision = false;
    uint32_t staleReason = SoBRLDatabaseSource::STALE_NONE;
};

template <typename Shape>
static void
apply_shape_material_policy(Shape &shape, const SoBRLDatabaseSource &source,
    const BObolMaterialPathState *resolved)
{
    if (resolved) {
	shape.materialColorValid = TRUE;
	shape.materialColor = resolved->color;
    } else if (source.databaseMetadataValid.getValue()) {
	shape.materialColorValid = source.databaseMaterialColorValid.getValue();
	shape.materialColor = source.databaseMaterialColor.getValue();
    }
    if (source.materialColorValid.getValue() &&
	(source.materialPolicy.getValue() != SoBRLDatabaseSource::MATERIAL_DATABASE ||
	 !shape.materialColorValid.getValue())) {
	shape.materialColorValid = TRUE;
	shape.materialColor = source.materialColor.getValue();
    }
    shape.materialRevision = source.materialRevision.getValue();
}

template <typename Shape>
static void
prepare_shape_material_policy(Shape &shape, const SoBRLDatabaseSource &source,
    BObolMaterialColorSweep &sweep)
{
    if (!shape.databaseIntent.getValue() && shape.nonDatabaseSource.getValue()) return;
    BObolMaterialPathState resolved;
    const bool found = sweep.resolve(shape.sourcePath.getValue().getString(), resolved);
    apply_shape_material_policy(shape, source, found ? &resolved : nullptr);
}

class PreparedOwnedMetadata {
public:
    template <typename Shape, typename Prepare>
    PreparedOwnedMetadata(Shape &shape, Prepare prepare) :
	node(&shape), candidate(new Shape), sharedGeometry(nullptr)
    {
	auto &next = *static_cast<Shape *>(this->candidate.get());
	copy_publication_scalar_fields(next, shape);
	prepare(next);
	this->prepareFields(shape);
    }

    template <typename Shape>
    PreparedOwnedMetadata(Shape &shape, const SoBRLDatabaseSource &source, bool clearGeometry) :
	PreparedOwnedMetadata(shape, source, source.sourceRevision.getValue(), clearGeometry)
    {}

    template <typename Shape>
    PreparedOwnedMetadata(Shape &shape, const SoBRLDatabaseSource &source,
	uint32_t revision, bool clearGeometry, OwnedShapeState state = OwnedShapeState::Realized,
	const char *oldSourcePath = nullptr, BObolMaterialColorSweep *materialSweep = nullptr) :
	node(&shape), candidate(new Shape),
	sharedGeometry(clearGeometry ? shape.sharedGeometry.getValue() : nullptr)
    {
	Shape &next = *static_cast<Shape *>(this->candidate.get());
	copy_publication_scalar_fields(next, shape);
	std::vector<SoField *> geometryFields;
	if (clearGeometry) {
	    const char *defaultType = std::is_same<Shape, SoBRLVListShape>::value ? "line-set" : "indexed-face-set";
	    const char *defaultKind = std::is_same<Shape, SoBRLVListShape>::value ? "line" : "surface";
	    assign_external_primary_identity(&next, &source,
		external_string_or_default(shape.sourceType.getValue().getString(), defaultType),
		external_string_or_default(shape.geometryKind.getValue().getString(), defaultKind));
	    clear_external_shape_geometry(next);
	    geometryFields = external_geometry_fields(shape);
	    this->clear = [](SoNode *target) { clear_external_shape_geometry(*static_cast<Shape *>(target)); };
	}
	if (state == OwnedShapeState::Realized)
	    prepare_realized_shape_owner(next, source, revision);
	else if (state == OwnedShapeState::Configured)
	    sync_shape_source_state(&next, &source);
	else if (state == OwnedShapeState::DisplayName)
	    sync_shape_display_name(&next, &source);
	else if (state == OwnedShapeState::Placement)
	    sync_shape_placement_state(&next, &source);
	else if (state == OwnedShapeState::Owner)
	    sync_shape_owner_state(&next, &source);
	else
	    prepare_shape_material_policy(next, source, *materialSweep);
	if (oldSourcePath)
	    retarget_realized_shape_source(&next, &source, oldSourcePath, source.path.getValue().getString(), revision);
	this->prepareFields(shape, geometryFields, clearGeometry ? &shape.sharedGeometry : nullptr);
    }
    PreparedOwnedMetadata(SoBRLMaterialObject &material, const char *oldSourcePath,
	const char *newSourcePath, uint32_t revision) :
	node(&material), candidate(new SoBRLMaterialObject), sharedGeometry(nullptr)
    {
	auto &next = *static_cast<SoBRLMaterialObject *>(this->candidate.get());
	copy_publication_scalar_fields(next, material);
	retarget_material_object_source(&next, oldSourcePath, newSourcePath, revision);
	this->prepareFields(material);
    }
    ~PreparedOwnedMetadata() { this->restore(); }
    bool changed() const { return this->notifications->changed(); }
    SoNode *target() const { return this->node.get(); }
    void commit()
    {
	if (this->clear) this->clear(this->node.get());
	this->values.commit();
    }
    void restore() { if (this->notifications) this->notifications->restore(); }
    void notify(std::exception_ptr &failure) { this->notifications->notify(failure); }
private:
    void prepareFields(SoNode &target, const std::vector<SoField *> &geometryFields = {},
	SoSFNode *clearedGeometry = nullptr)
    {
	const SoFieldData *data = static_cast<const SoFieldContainer &>(target).getFieldData();
	std::vector<PublicationFieldChange> changed;
	changed.reserve(size_t(data->getNumFields()));
	this->values.reserve(size_t(data->getNumFields()));
	for (int i = 0; i < data->getNumFields(); ++i) {
	    SoField *field = data->getField(&target, i);
	    const SoField *value = this->candidate->getField(data->getFieldName(i));
	    bool differs = this->values.prepare(*field, value);
	    if (!differs && std::find(geometryFields.begin(), geometryFields.end(), field) != geometryFields.end())
		differs = static_cast<SoMField *>(field)->getNum() != 0;
	    else if (!differs && field == clearedGeometry)
		differs = clearedGeometry->getValue() != nullptr;
	    changed.push_back({field, differs});
	}
	this->notifications = std::make_unique<PreparedNotifications<std::vector<PublicationFieldChange>>>(
	    target, std::move(changed));
    }
    SbModernUtils::SoNodeRef node;
    SbModernUtils::SoNodeRef candidate;
    SbModernUtils::SoNodeRef sharedGeometry;
    PreparedScalarValues values;
    std::unique_ptr<PreparedNotifications<std::vector<PublicationFieldChange>>> notifications;
    void (*clear)(SoNode *) = nullptr;
};

template <typename Visit>
static void
visit_owned_nodes(SoNode *node, std::unordered_set<SoNode *> &visited, Visit visit)
{
    if (!node || node->isOfType(SoBRLDatabaseSource::getClassTypeId()) || !visited.insert(node).second)
	return;
    visit(node);
    if (node->isOfType(SoGroup::getClassTypeId())) {
	auto *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); ++i) visit_owned_nodes(group->getChild(i), visited, visit);
    }
}

static void
prepare_owned_shape_edit(SoNode *node, const SoBRLDatabaseSource &source,
    uint32_t revision, std::vector<std::unique_ptr<PreparedOwnedMetadata>> &shapes,
    std::unordered_set<SoNode *> &visited, OwnedShapeState state = OwnedShapeState::Realized,
    const char *oldSourcePath = nullptr, BObolMaterialColorSweep *materialSweep = nullptr)
{
    visit_owned_nodes(node, visited, [&](SoNode *owned) {
	if (owned->isOfType(SoBRLVListShape::getClassTypeId()))
	    shapes.push_back(std::make_unique<PreparedOwnedMetadata>(
		*static_cast<SoBRLVListShape *>(owned), source, revision, false, state, oldSourcePath, materialSweep));
	else if (owned->isOfType(SoBRLMeshShape::getClassTypeId()))
	    shapes.push_back(std::make_unique<PreparedOwnedMetadata>(
		*static_cast<SoBRLMeshShape *>(owned), source, revision, false, state, oldSourcePath, materialSweep));
	else if (oldSourcePath && owned->isOfType(SoBRLMaterialObject::getClassTypeId()))
	    shapes.push_back(std::make_unique<PreparedOwnedMetadata>(
		*static_cast<SoBRLMaterialObject *>(owned), oldSourcePath, source.path.getValue().getString(), revision));
    });
}

struct BObolSourceRenamePublication::Impl {
    Impl(SoBRLDatabaseSource &target, const char *oldObjectPath,
	const char *newObjectPath, const char *newSourceInstanceKey,
	const char *newSourcePath, uint32_t requestedRevision,
	const char *newParentInstanceKey) :
	source(target), retainedSource(&target), owner(new SoBRLDatabaseSource(
	    SoBRLDatabaseSource::FieldObservation::Detached)),
	next(*static_cast<SoBRLDatabaseSource *>(this->owner.get()))
    {
	this->next.copyFieldValues(&target, FALSE);

	if (newSourcePath && newSourcePath[0]) {
	    const char *nextInstance = newSourceInstanceKey &&
		newSourceInstanceKey[0] ? newSourceInstanceKey : newSourcePath;
	    const bool pathChanged = !database_source_string_equal(
		target.path.getValue(), newSourcePath);
	    const bool instanceChanged = !database_source_string_equal(
		target.instanceKey.getValue(), nextInstance);
	    uint32_t revision = requestedRevision;
	    if (!revision) {
		revision = target.sourceRevision.getValue();
		if (pathChanged || instanceChanged)
		    revision = bobol_identity_successor_or_terminate(revision);
	    }
	    const bool revisionChanged =
		target.sourceRevision.getValue() != revision;
	    this->identityChanged = pathChanged || instanceChanged ||
		revisionChanged;
	    if (this->identityChanged) {
		const SbString oldInstance = target.instanceKey.getValue();
		const SbString oldPath = target.path.getValue();
		this->next.instanceKey = nextInstance;
		this->next.path = newSourcePath;
		this->next.sourceRevision = revision;
		if (database_source_string_equal(
			target.representationKey.getValue(),
			oldInstance.getString()))
		    this->next.representationKey = nextInstance;
		else if (database_source_string_equal(
			target.representationKey.getValue(),
			oldPath.getString()))
		    this->next.representationKey = newSourcePath;
		this->next.markStale(SoBRLDatabaseSource::STALE_SOURCE);
		this->next.realizationIdentity =
		    source_realization_identity(&this->next);

		std::unordered_set<SoNode *> visited;
		for (int i = 0; i < target.getNumChildren(); ++i)
		    prepare_owned_shape_edit(target.getChild(i), this->next,
			revision, this->shapes, visited,
			OwnedShapeState::Configured, oldPath.getString());
	    }
	}
	if (newParentInstanceKey && !database_source_string_equal(
		target.parentInstanceKey.getValue(), newParentInstanceKey)) {
	    this->next.parentInstanceKey = newParentInstanceKey;
	    this->hierarchyChanged = true;
	}
	if (this->identityChanged || this->hierarchyChanged) {
	    this->fields = std::make_unique<PreparedScalarFields>(target,
		this->next);
	}

	if (target.d->compactIndex) {
	    auto candidate = std::make_unique<BObolCompactInstanceIndex>(
		*target.d->compactIndex);
	    this->compactChanges = compact_retarget_occurrence_paths(*candidate,
		oldObjectPath, newObjectPath);
	    if (this->hierarchyChanged) {
		const Obol::InstanceId parent = cad_source_parent_instance(
		    &this->next);
		for (auto &instance : candidate->instances)
		    instance.record.parent = parent;
	    }
	    if (this->compactChanges > 0 || this->hierarchyChanged) {
		this->compact = std::move(candidate);
		this->compactChanged = true;
	    }
	}

	if (this->compactChanged && !this->fields) {
	    this->nodeNotifications = std::make_unique<PreparedNotifications<
		std::vector<PublicationFieldChange>>>(target,
		std::vector<PublicationFieldChange>());
	}
    }

    bool changed() const
    {
	return this->identityChanged || this->hierarchyChanged ||
	    this->compactChanged;
    }
    int compactPathChangeCount() const { return this->compactChanges; }

    void commit() noexcept
    {
	if (!this->changed())
	    return;
	if (this->fields)
	    this->fields->commit();
	for (auto &shape : this->shapes)
	    shape->commit();
	if (this->compact) {
	    this->previousCompact.reset(this->source.d->compactIndex);
	    this->source.d->compactIndex = this->compact.release();
	}
	if (this->identityChanged)
	    this->source.revokeSourceResourceAuthority();
	if (this->hierarchyChanged)
	    this->source.d->clearCompiledCompactEvidence();
	this->source.markCompiledAssemblyDirty();
	this->source.markCadBatchDirty();
	this->source.markDisplayMeshLodDirty();
    }

    void restore()
    {
	if (this->fields)
	    this->fields->restore();
	if (this->nodeNotifications)
	    this->nodeNotifications->restore();
	for (auto &shape : this->shapes)
	    shape->restore();
    }

    void notify(std::exception_ptr &failure)
    {
	this->restore();
	if (this->fields)
	    this->fields->notify(failure);
	else if (this->nodeNotifications)
	    this->nodeNotifications->notify(failure);
	for (auto &shape : this->shapes)
	    shape->notify(failure);
    }

    SoBRLDatabaseSource &source;
    SbModernUtils::SoNodeRef retainedSource;
    SbModernUtils::SoNodeRef owner;
    SoBRLDatabaseSource &next;
    std::unique_ptr<PreparedScalarFields> fields;
    std::unique_ptr<PreparedNotifications<
	std::vector<PublicationFieldChange>>> nodeNotifications;
    std::vector<std::unique_ptr<PreparedOwnedMetadata>> shapes;
    std::unique_ptr<BObolCompactInstanceIndex> compact;
    std::unique_ptr<BObolCompactInstanceIndex> previousCompact;
    bool identityChanged = false;
    bool hierarchyChanged = false;
    bool compactChanged = false;
    int compactChanges = 0;
};

BObolSourceRenamePublication::BObolSourceRenamePublication(
    SoBRLDatabaseSource &source, const char *oldObjectPath,
    const char *newObjectPath, const char *newSourceInstanceKey,
    const char *newSourcePath, uint32_t sourceRevision,
    const char *newParentInstanceKey) :
    d(std::make_unique<Impl>(source, oldObjectPath, newObjectPath,
	newSourceInstanceKey, newSourcePath, sourceRevision,
	newParentInstanceKey))
{
}

BObolSourceRenamePublication::~BObolSourceRenamePublication() = default;

bool
BObolSourceRenamePublication::changed(void) const
{
    return this->d->changed();
}

int
BObolSourceRenamePublication::compactPathChangeCount(void) const
{
    return this->d->compactPathChangeCount();
}

void
BObolSourceRenamePublication::commit(void) noexcept
{
    this->d->commit();
}

void
BObolSourceRenamePublication::restore(void)
{
    this->d->restore();
}

void
BObolSourceRenamePublication::notify(std::exception_ptr &failure)
{
    this->d->notify(failure);
}

/* Database refresh changes intrinsic metadata and RGB only. Selection,
 * visibility, line styling and presentation opacity have separate owners. */
class PreparedCompactMaterials {
public:
    PreparedCompactMaterials(BObolCompactInstanceIndex *target, BObolMaterialColorSweep &sweep) : index(target)
    {
	if (!target) return;
	for (size_t i = 0; i < target->entries.size(); ++i) {
	    const auto &semantic = target->entries[i].semantic;
	    BObolMaterialPathState resolved;
	    if (!sweep.resolve(semantic.path.getString(), resolved)) continue;
	    if (semantic.regionId == resolved.regionId && semantic.airCode == resolved.airCode &&
		semantic.materialId == resolved.materialId && semantic.los == resolved.los &&
		semantic.materialColorValid && database_source_color_equal(semantic.materialColor, resolved.color) &&
		resolved.shader == semantic.materialShader.getString()) continue;
	    this->changes.push_back({i, BObolCompactMetadataEdit(semantic, resolved.regionId, resolved.airCode,
		resolved.materialId, resolved.los, TRUE, resolved.color, resolved.shader.c_str())});
	}
    }
    bool changed() const { return !this->changes.empty(); }
    void commit(const SoBRLDatabaseSource &source)
    {
	if (!this->index) return;
	for (auto &change : this->changes) change.metadata.commit(this->index->entries[change.ordinal]);
	size_t changedOrdinal = 0;
	for (size_t i = 0; i < this->index->entries.size(); ++i) {
	    auto &entry = this->index->entries[i];
	    const bool semanticChanged = changedOrdinal < this->changes.size() && this->changes[changedOrdinal].ordinal == i;
	    if (semanticChanged) ++changedOrdinal;
	    if (compact_set_material_styles(entry, source, true) || semanticChanged)
		entry.appearanceRevision = compact_next_revision(entry.appearanceRevision);
	    entry.shapeSummary.materialRevision = source.materialRevision.getValue();
	    entry.style = compact_effective_style(entry);
	    compact_sync_instance_style(*this->index, i);
	}
    }
private:
    struct Change {
	size_t ordinal;
	BObolCompactMetadataEdit metadata;
    };
    BObolCompactInstanceIndex *index;
    std::vector<Change> changes;
};

template <typename Prepare>
static std::vector<std::unique_ptr<PreparedOwnedMetadata>>
prepare_owned_metadata(SoBRLDatabaseSource &source, Prepare prepare)
{
    std::vector<std::unique_ptr<PreparedOwnedMetadata>> shapes;
    std::unordered_set<SoNode *> visited;
    for (int i = 0; i < source.getNumChildren(); ++i) {
	visit_owned_nodes(source.getChild(i), visited, [&](SoNode *owned) {
	    std::unique_ptr<PreparedOwnedMetadata> edit;
	    if (owned->isOfType(SoBRLVListShape::getClassTypeId()))
		edit = std::make_unique<PreparedOwnedMetadata>(*static_cast<SoBRLVListShape *>(owned), prepare);
	    else if (owned->isOfType(SoBRLMeshShape::getClassTypeId()))
		edit = std::make_unique<PreparedOwnedMetadata>(*static_cast<SoBRLMeshShape *>(owned), prepare);
	    if (edit && edit->changed()) shapes.push_back(std::move(edit));
	});
    }
    return shapes;
}

static bool
compact_presentation_rule_equal(
    const BObolCompactOccurrenceRegistryState::PresentationOverride &left,
    const BObolCompactOccurrenceRegistryState::PresentationOverride &right)
{
    using Rule = BObolCompactOccurrenceRegistryState::PresentationOverride;
    return compact_presentation_override_same_key(left, right) &&
	(left.property == Rule::TRANSPARENCY ?
	 !database_source_float_different(left.transparency,
	     right.transparency) : left.state == right.state);
}

static void
compact_presentation_rule_set(
    std::vector<BObolCompactOccurrenceRegistryState::PresentationOverride> &rules,
    const BObolCompactPresentationPatch &patch,
    BObolCompactOccurrenceRegistryState::PresentationOverride::Property property)
{
    using Rule = BObolCompactOccurrenceRegistryState::PresentationOverride;
    Rule next;
    next.property = property;
    next.path = database_source_skip_leading_slash(patch.path.getString());
    next.match = patch.match;
    next.state = property == Rule::VISIBILITY ? patch.visible :
	patch.highlighted;
    next.transparency = std::max(0.0f,
	std::min(1.0f, patch.transparency));
    if (!rules.empty() && compact_presentation_rule_equal(rules.back(), next))
	return;
    rules.erase(std::remove_if(rules.begin(), rules.end(),
	[&next](const Rule &rule) {
	    return compact_presentation_override_same_key(rule, next);
	}), rules.end());
    rules.push_back(std::move(next));
}

static bool
compact_presentation_rules_equal(
    const std::vector<BObolCompactOccurrenceRegistryState::PresentationOverride> &left,
    const std::vector<BObolCompactOccurrenceRegistryState::PresentationOverride> &right)
{
    return left.size() == right.size() &&
	std::equal(left.begin(), left.end(), right.begin(),
	    compact_presentation_rule_equal);
}

SbBool
SoBRLDatabaseSource::presentationPatchValid(
    const BObolDatabaseSourcePresentationPatch &patch)
{
    for (const auto &compact : patch.compactOverrides) {
	if (!compact.visibleValid && !compact.highlightedValid &&
	    !compact.transparencyValid)
	    continue;
	if (compact.match != BOBOL_COMPACT_PATH_EXACT &&
	    compact.match != BOBOL_COMPACT_PATH_SUBTREE &&
	    compact.match != BOBOL_COMPACT_PATH_OBJECT)
	    return false;
    }
    return true;
}

/* Only membership buffers, retained-rule strings and a frontier mask need
 * allocation. Style and selection changes commit in place without copying
 * geometry or the registry. */
class PreparedCompactDisplay {
public:
    PreparedCompactDisplay(BObolCompactInstanceIndex *target, const SoBRLDatabaseSource &next,
	const SourceDisplayChanges &edit, BObolCompactOccurrenceRegistryState &presentationState,
	bool registry) : index(target), source(next), changes(edit), presentation(presentationState),
	select(!registry && edit.patch.selectedValid),
	highlight(!registry && edit.patch.highlightedValid)
    {
	const auto *intent = edit.presentation;
	if (intent) {
	    this->rules = presentationState.compactPresentationOverrides;
	    if (intent->clearCompactHighlightOverrides) {
		using Rule = BObolCompactOccurrenceRegistryState::PresentationOverride;
		this->rules.erase(std::remove_if(this->rules.begin(), this->rules.end(),
		    [](const Rule &rule) { return rule.property == Rule::HIGHLIGHT; }),
		    this->rules.end());
	    }
	    for (const auto &patch : intent->compactOverrides) {
		if (patch.visibleValid) {
		    compact_presentation_rule_set(this->rules, patch,
			BObolCompactOccurrenceRegistryState::PresentationOverride::VISIBILITY);
		    this->visibilityRulesEdited = true;
		}
		if (patch.highlightedValid)
		    compact_presentation_rule_set(this->rules, patch,
			BObolCompactOccurrenceRegistryState::PresentationOverride::HIGHLIGHT);
		if (patch.transparencyValid)
		    compact_presentation_rule_set(this->rules, patch,
			BObolCompactOccurrenceRegistryState::PresentationOverride::TRANSPARENCY);
	    }
	    this->rulesChanged = !compact_presentation_rules_equal(this->rules,
		presentationState.compactPresentationOverrides);
	}
	if (!this->index) {
	    this->compactChanged = this->rulesChanged;
	    return;
	}
	this->evaluateVisibility = edit.patch.visibleValid ||
	    this->visibilityRulesEdited;
	if (this->evaluateVisibility) {
	    if (presentationState.compactVisibilityFrontierActive)
		this->allowed = compact_visibility_frontier_mask(*target, presentationState, 0);
	    this->hidden.reserve(target->entries.size());
	    this->visibilityChanges.reserve(target->entries.size());
	    for (size_t i = 0; i < target->entries.size(); ++i) {
		const auto &entry = target->entries[i];
		const BObolCompactPresentationOverrideState state =
		    this->nextPresentation(entry);
		const SbBool visible = this->nextVisible(entry, state, i);
		if (!visible) this->hidden.push_back(entry.instance);
		if (visible != entry.visible) {
		    this->visibilityChanges.push_back(i);
		    this->compactChanged = true;
		}
	    }
	}
	if (this->select && edit.patch.selected) {
	    this->selected.reserve(target->entries.size());
	    for (const auto &entry : target->entries) this->selected.push_back(entry.instance);
	}
	this->compactChanged = this->compactChanged || this->rulesChanged;
	if (intent && intent->allCompactHighlightedValid) {
	    for (const auto &entry : target->entries)
		this->compactChanged = this->compactChanged ||
		    entry.authoredHighlighted != intent->allCompactHighlighted;
	}
    }
    void commit()
    {
	if (this->rulesChanged) {
	    auto &current = this->presentation.compactPresentationOverrides;
	    current.swap(this->rules);
	}
	if (!this->index) return;
	const auto &patch = this->changes.patch;
	for (size_t i = 0; i < this->index->entries.size(); ++i) {
	    auto &entry = this->index->entries[i];
	    const SbBool oldTransparencyValid = entry.presentationTransparencyValid;
	    const float oldTransparency = entry.presentationTransparency;
	    if (this->rulesChanged) {
		const BObolCompactPresentationOverrideState state =
		    compact_presentation_override_state(entry,
			this->presentation.compactPresentationOverrides);
		entry.presentationVisibleValid = state.visibleValid;
		entry.presentationVisible = state.visible;
		entry.presentationHighlightedValid = state.highlightedValid;
		entry.presentationHighlighted = state.highlighted;
		entry.presentationTransparencyValid = state.transparencyValid;
		entry.presentationTransparency = state.transparency;
	    }
	    if (patch.visibleValid) {
		entry.authoredVisible = this->nextAuthoredVisible(entry);
	    }
	    if (this->evaluateVisibility) {
		BObolCompactPresentationOverrideState state;
		state.visibleValid = entry.presentationVisibleValid;
		state.visible = entry.presentationVisible;
		const SbBool visible = this->nextVisible(entry, state, i);
		if (entry.visible != visible) entry.visibilityRevision = compact_next_revision(entry.visibilityRevision);
		entry.visible = visible;
	    }
	    const SbBool previousSelected = entry.selected;
	    const SbBool previousHighlighted = entry.highlighted;
	    if (this->select) entry.selected = patch.selected;
	    if (this->highlight)
		entry.authoredHighlighted = patch.highlighted;
	    if (this->changes.presentation &&
		this->changes.presentation->allCompactHighlightedValid)
		entry.authoredHighlighted =
		    this->changes.presentation->allCompactHighlighted;
	    entry.highlighted = compact_effective_highlight(entry);
	    if (entry.selected != previousSelected || entry.highlighted != previousHighlighted)
		entry.selectionRevision = compact_next_revision(entry.selectionRevision);
	    bool appearanceChanged = this->updateStyle(entry.normalStyle, entry, FALSE, FALSE);
	    appearanceChanged = this->updateStyle(entry.selectedStyle, entry, TRUE, FALSE) || appearanceChanged;
	    appearanceChanged = this->updateStyle(entry.highlightedStyle, entry, FALSE, TRUE) || appearanceChanged;
	    const Obol::InstanceStyle effectiveStyle = compact_effective_style(entry);
	    appearanceChanged = !compact_style_equal(entry.style, effectiveStyle) ||
		oldTransparencyValid != entry.presentationTransparencyValid ||
		(oldTransparencyValid && entry.presentationTransparencyValid &&
		 database_source_float_different(oldTransparency,
		     entry.presentationTransparency)) || appearanceChanged;
	    if (appearanceChanged) entry.appearanceRevision = compact_next_revision(entry.appearanceRevision);
	    entry.style = effectiveStyle;
	    compact_sync_instance_style(*this->index, i);
	    compact_sync_shape_display_summary(entry);
	    if (patch.colorOverrideValid) entry.shapeSummary.colorOverride = patch.colorOverride;
	    if (patch.colorValid) entry.shapeSummary.color = patch.color;
	    if (this->changes.materialRevision) entry.shapeSummary.materialRevision = this->source.materialRevision.getValue();
	}
	if (this->evaluateVisibility) this->index->hiddenInstances.swap(this->hidden);
	if (this->select) this->index->selectedInstances.swap(this->selected);
    }
    bool changed() const { return this->compactChanged; }
    const std::vector<size_t> &visibilityChangedEntries() const { return this->visibilityChanges; }
private:
    SbBool nextAuthoredVisible(const BObolCompactInstanceEntry &entry) const
    {
	// Retired whole-target overviews cannot be resurrected by a display edit.
	if (compact_retired_overview(entry)) return FALSE;
	return this->changes.patch.visibleValid ? this->changes.patch.visible :
	    entry.authoredVisible;
    }
    BObolCompactPresentationOverrideState nextPresentation(
	const BObolCompactInstanceEntry &entry) const
    {
	return compact_presentation_override_state(entry,
	    this->changes.presentation ? this->rules :
	    this->presentation.compactPresentationOverrides);
    }
    SbBool nextVisible(const BObolCompactInstanceEntry &entry,
	const BObolCompactPresentationOverrideState &state,
	size_t ordinal) const
    {
	if (compact_retired_overview(entry)) return FALSE;
	return (state.visibleValid ? state.visible : this->nextAuthoredVisible(entry)) &&
	    (this->allowed.empty() || this->allowed[ordinal]);
    }
    bool updateStyle(Obol::InstanceStyle &style, const BObolCompactInstanceEntry &entry,
	SbBool selectedState, SbBool highlightedState) const
    {
	const auto &patch = this->changes.patch;
	bool colorChanged = patch.colorOverrideValid || patch.colorValid || this->changes.materialColor;
	if (highlightedState) colorChanged = patch.highlightedColorValid;
	else if (selectedState) colorChanged = patch.selectedColorValid;
	else if (entry.shapeSummary.ghosted) colorChanged = patch.ghostedColorValid;
	if (!colorChanged && !patch.transparencyValid && !patch.lineWidthValid && !patch.lineStyleValid) return false;
	const auto requested = compact_entry_style_from_source(&this->source, entry, selectedState, highlightedState);
	auto next = style;
	if (colorChanged) {
	    next.hasColorOverride = requested.hasColorOverride;
	    for (int axis = 0; axis < 3; ++axis) next.color[axis] = requested.color[axis];
	}
	if (patch.transparencyValid) next.color[3] = requested.color[3];
	if (patch.lineWidthValid) next.lineWidth = requested.lineWidth;
	if (patch.lineStyleValid) next.linePattern = requested.linePattern;
	if (compact_style_equal(style, next)) return false;
	style = next;
	return true;
    }
    BObolCompactInstanceIndex *index;
    const SoBRLDatabaseSource &source;
    const SourceDisplayChanges &changes;
    BObolCompactOccurrenceRegistryState &presentation;
    bool select, highlight;
    bool rulesChanged = false;
    bool visibilityRulesEdited = false;
    bool evaluateVisibility = false;
    bool compactChanged = false;
    std::vector<BObolCompactOccurrenceRegistryState::PresentationOverride> rules;
    std::vector<SbBool> allowed;
    std::vector<Obol::InstanceId> hidden, selected;
    std::vector<size_t> visibilityChanges;
};

static void
compact_sync_draw_mode(BObolCompactInstanceIndex *index, const SoBRLDatabaseSource &source)
{
    if (!index) return;
    const int mode = source_record_draw_mode(&source);
    for (auto &entry : index->entries) {
	entry.shapeSummary.drawMode = mode;
	entry.shapeSummary.hiddenLine = mode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    }
}

static void
compact_sync_material_policy(BObolCompactInstanceIndex *index, const SoBRLDatabaseSource &source)
{
    if (!index) return;
    for (size_t i = 0; i < index->entries.size(); ++i) {
	auto &entry = index->entries[i];
	if (compact_set_material_styles(entry, source, true))
	    entry.appearanceRevision = compact_next_revision(entry.appearanceRevision);
	entry.style = compact_effective_style(entry);
	compact_sync_instance_style(*index, i);
    }
}

class SoBRLDatabaseSource::ConfigurationPublication {
public:
    explicit ConfigurationPublication(SoBRLDatabaseSource &target) :
	source(target), owner(new SoBRLDatabaseSource(FieldObservation::Detached)),
	next(*static_cast<SoBRLDatabaseSource *>(this->owner.get())), database(target.d->dbip)
    {
	this->next.copyFieldValues(&target, FALSE);
    }
    void publish(uint32_t reason, bool retargetIdentity = false,
	PublicationCommit committed = nullptr, void *context = nullptr)
    {
	this->next.markStale(reason);
	if (retargetIdentity)
	    this->next.realizationIdentity = source_realization_identity(&this->next);
	std::vector<std::unique_ptr<PreparedOwnedMetadata>> shapes;
	std::unordered_set<SoNode *> visited;
	const char *oldPath = retargetIdentity ? this->source.path.getValue().getString() : nullptr;
	for (int i = 0; i < this->source.getNumChildren(); ++i)
	    prepare_owned_shape_edit(this->source.getChild(i), this->next, this->next.sourceRevision.getValue(),
		shapes, visited, OwnedShapeState::Configured, oldPath);
	PreparedScalarFields fields(this->source, this->next);

	(void)this->source.d->setDatabaseBinding(this->database);
	fields.commit();
	for (auto &shape : shapes) shape->commit();
	if (invalidation_revokes_source_contract(reason))
	    this->source.revokeSourceResourceAuthority();
	// Configuration changes draw channels, not occurrence identity or style.
	compact_sync_draw_mode(this->source.d->compactIndex, this->source);
	this->finish(fields, shapes, committed, context);
    }
    bool publishDisplay(const SourceDisplayChanges &changes, PublicationCommit committed, void *context)
    {
	this->next.markStale(changes.staleReason);
	BObolMaterialColorSweep sweep(this->database);
	const auto prepare = [&](auto &shape) {
	    if (changes.staleReason) sync_shape_owner_state(&shape, &this->next);
	    if (!shape.databaseIntent.getValue() && shape.nonDatabaseSource.getValue()) return;
	    apply_display_patch(shape, changes.patch);
	    if (changes.materialColor) prepare_shape_material_policy(shape, this->next, sweep);
	    else if (changes.materialRevision) shape.materialRevision = this->next.materialRevision.getValue();
	};
	auto shapes = prepare_owned_metadata(this->source, prepare);
	PreparedCompactDisplay compact(this->source.d->compactIndex, this->next, changes,
	    *this->source.d, this->source.isCompactOccurrenceRegistry());
	PreparedScalarFields fields(this->source, this->next);
	if (!fields.changed() && shapes.empty() && !compact.changed())
	    return false;
	fields.commit();
	for (auto &shape : shapes) shape->commit();
	compact.commit();
	if (invalidation_revokes_source_contract(changes.staleReason))
	    this->source.revokeSourceResourceAuthority();
	this->source.markDisplayMeshLodVisibilityDirty(compact.visibilityChangedEntries());
	this->finish(fields, shapes, committed, context);
	return true;
    }
    void publishDatabaseMetadata(PublicationCommit committed, void *context)
    {
	const std::string sourcePath = database_source_db_path_without_instance_suffixes(this->next.path.getValue().getString());
	auto shapes = prepare_owned_metadata(this->source, [&](auto &shape) {
	    if (!shape.databaseIntent.getValue() && shape.nonDatabaseSource.getValue()) return;
	    // Aggregate metadata cannot describe a different descendant's region or shader.
	    if (shape.sourcePath.getValue().getLength() &&
		database_source_db_path_without_instance_suffixes(shape.sourcePath.getValue().getString()) != sourcePath) return;
	    apply_source_database_metadata(&shape, &this->next);
	});
	PreparedScalarFields fields(this->source, this->next);
	fields.commit();
	for (auto &shape : shapes) shape->commit();
	this->finish(fields, shapes, committed, context);
    }
    int refreshMaterials(BObolMaterialColorSweep &sweep, uint32_t revision,
	PublicationCommit committed, void *context)
    {
	const bool compactSource = this->source.d->compactIndex != nullptr;
	BObolMaterialPathState resolved;
	const bool sourceResolved = sweep.resolve(this->next.path.getValue().getString(), resolved);
	if (!sourceResolved && !compactSource) return 0;
	bool changed = false;
	if (sourceResolved) {
	    changed = !this->next.materialColorValid.getValue() ||
		!database_source_color_equal(this->next.materialColor.getValue(), resolved.color);
	    this->next.materialColorValid = TRUE;
	    if (changed) this->next.materialColor = resolved.color;
	    if (!compactSource && resolved.inRegion) {
		this->next.databaseMetadataValid = TRUE;
		this->next.databaseRegionId = resolved.regionId;
		this->next.databaseAirCode = resolved.airCode;
		this->next.databaseMaterialId = resolved.materialId;
		this->next.databaseLos = resolved.los;
		this->next.databaseMaterialColorValid = TRUE;
		this->next.databaseMaterialColor = resolved.color;
		this->next.databaseMaterialShader = resolved.shader.c_str();
	    }
	}
	PreparedCompactMaterials compact(this->source.d->compactIndex, sweep);
	changed = changed || compact.changed();
	this->next.materialRevision = revision;
	auto shapes = prepare_owned_metadata(this->source, [&](auto &shape) {
	    if (!shape.databaseIntent.getValue() && shape.nonDatabaseSource.getValue()) return;
	    BObolMaterialPathState material;
	    const SbString &shapePath = shape.sourcePath.getValue().getLength() ? shape.sourcePath.getValue() : this->next.path.getValue();
	    const bool found = sweep.resolve(shapePath.getString(), material);
	    if (found) {
		changed = changed || shape.regionId.getValue() != material.regionId || shape.airCode.getValue() != material.airCode ||
		    shape.materialId.getValue() != material.materialId || shape.los.getValue() != material.los ||
		    material.shader != shape.materialShader.getValue().getString();
		shape.regionId = material.regionId;
		shape.airCode = material.airCode;
		shape.materialId = material.materialId;
		shape.los = material.los;
		shape.materialShader = material.shader.c_str();
	    }
	    const SbBool colorValid = shape.materialColorValid.getValue();
	    const SbColor oldColor = shape.materialColor.getValue();
	    apply_shape_material_policy(shape, this->next, found ? &material : nullptr);
	    changed = changed || colorValid != shape.materialColorValid.getValue() ||
		!database_source_color_equal(oldColor, shape.materialColor.getValue());
	});
	// Compact stamp-only refreshes acknowledge the sweep without invalidating presentation.
	if (compactSource && !changed) {
	    shapes.clear();
	    PreparedFieldNotifications<1> fields(this->source, {{{&this->source.materialRevision, true}}});
	    this->source.materialRevision = revision;
	    std::exception_ptr failure;
	    fields.notify(failure);
	    if (failure) std::rethrow_exception(failure);
	    return 0;
	}
	PreparedScalarFields fields(this->source, this->next);
	fields.commit();
	for (auto &shape : shapes) shape->commit();
	compact.commit(this->next);
	this->finish(fields, shapes, committed, context);
	return 1;
    }
    void publishViewPolicy(PublicationCommit committed, void *context)
    {
	this->next.markStale(STALE_VIEW);
	std::vector<std::unique_ptr<PreparedOwnedMetadata>> shapes;
	std::unordered_set<SoNode *> visited;
	for (int i = 0; i < this->source.getNumChildren(); ++i)
	    prepare_owned_shape_edit(this->source.getChild(i), this->next, this->next.sourceRevision.getValue(),
		shapes, visited, OwnedShapeState::Owner);
	PreparedScalarFields fields(this->source, this->next);
	fields.commit();
	for (auto &shape : shapes) shape->commit();
	this->finish(fields, shapes, committed, context);
    }
    void publishMetadata(MetadataIntent intent, PublicationCommit committed, void *context)
    {
	std::vector<std::unique_ptr<PreparedOwnedMetadata>> shapes;
	if (intent != MetadataIntent::Hierarchy) {
	    std::unordered_set<SoNode *> visited;
	    BObolMaterialColorSweep sweep(this->database);
	    const OwnedShapeState state = intent == MetadataIntent::DisplayName ?
		OwnedShapeState::DisplayName : OwnedShapeState::MaterialPolicy;
	    for (int i = 0; i < this->source.getNumChildren(); ++i)
		prepare_owned_shape_edit(this->source.getChild(i), this->next, this->next.sourceRevision.getValue(),
		    shapes, visited, state, nullptr, &sweep);
	    shapes.erase(std::remove_if(shapes.begin(), shapes.end(),
		[](const std::unique_ptr<PreparedOwnedMetadata> &shape) { return !shape->changed(); }), shapes.end());
	}
	std::vector<std::pair<SbString *, SbString>> names;
	auto *index = this->source.d->compactIndex;
	if (index && intent == MetadataIntent::DisplayName) {
	    for (auto &entry : index->entries) {
		auto &summary = entry.shapeSummary;
		if (!summary.databaseIntent && summary.nonDatabaseSource) continue;
		const SbString &name = this->next.displayName.getValue().getLength() ?
		    this->next.displayName.getValue() : (summary.sourceName.getLength() ? summary.sourceName : summary.path);
		if (summary.displayName != name) names.emplace_back(&summary.displayName, name);
	    }
	}
	const Obol::InstanceId parent = index && intent == MetadataIntent::Hierarchy ?
	    cad_source_parent_instance(&this->next) : Obol::CadIdBuilder::rootInstance();
	PreparedScalarFields fields(this->source, this->next);
	fields.commit();
	for (auto &shape : shapes) shape->commit();
	for (auto &name : names) *name.first = std::move(name.second);
	if (intent == MetadataIntent::MaterialPolicy) compact_sync_material_policy(index, this->source);
	if (index && intent == MetadataIntent::Hierarchy) {
	    for (auto &instance : index->instances) instance.record.parent = parent;
	    // Parent identity is outside the per-leaf geometry/style stamps.
	    this->source.d->clearCompiledCompactEvidence();
	}
	this->finish(fields, shapes, committed, context);
    }
private:
    void finish(PreparedScalarFields &fields, std::vector<std::unique_ptr<PreparedOwnedMetadata>> &shapes,
	PublicationCommit committed, void *context)
    {
	this->source.markCompiledAssemblyDirty();
	this->source.markCadBatchDirty();
	if (committed) committed(context);
	fields.restore();
	for (auto &shape : shapes) shape->restore();
	std::exception_ptr failure;
	fields.notify(failure);
	for (auto &shape : shapes) shape->notify(failure);
	if (failure) std::rethrow_exception(failure);
    }
    SoBRLDatabaseSource &source;
    SbModernUtils::SoNodeRef owner;
public:
    SoBRLDatabaseSource &next;
    struct db_i *database;
};

void
SoBRLDatabaseSource::failSafeObservedFieldChange(uint32_t reason) noexcept
{
    /* The public field has committed before notify() can prepare its dependent
     * publication.  Establish an allocation-free stale state first so a later
     * preparation failure cannot leave changed identity authorizing old
     * resources.  A repeated field notification completes the owned records. */
    const bool retainExternalFailure = reason == STALE_VIEW &&
	(this->realizationRoleFlags.getValue() & REALIZATION_ROLE_EXTERNAL) &&
	this->realizationStatus.getValue() == FAILED;
    const SbBool statusNotifications = this->realizationStatus.enableNotify(FALSE);
    const SbBool diagnosticNotifications = this->realizationDiagnostic.enableNotify(FALSE);
    const SbBool staleNotifications = this->stale.enableNotify(FALSE);
    const SbBool reasonNotifications = this->staleReason.enableNotify(FALSE);
    const SbBool boundsNotifications = this->sourceBoundsExact.enableNotify(FALSE);

    if (!retainExternalFailure) {
	this->realizationStatus = UNREALIZED;
	this->realizationDiagnostic.setValue(SbString());
    }
    this->stale = TRUE;
    this->staleReason = this->staleReason.getValue() | reason;
    if (invalidation_revokes_source_contract(reason)) {
	this->revokeSourceResourceAuthority();
	this->sourceBoundsExact = FALSE;
    }
    this->markCompiledAssemblyDirty();
    this->markCadBatchDirty();

    this->realizationStatus.enableNotify(statusNotifications);
    this->realizationDiagnostic.enableNotify(diagnosticNotifications);
    this->stale.enableNotify(staleNotifications);
    this->staleReason.enableNotify(reasonNotifications);
    this->sourceBoundsExact.enableNotify(boundsNotifications);
}

void
SoBRLDatabaseSource::publishObservedFieldChange(uint32_t reason)
{
    /* The triggering field has already changed. Complete every dependent
     * source and owned-shape record before that notification continues to
     * later auditors. Reusing the configuration publication also keeps
     * invalidation and resource retirement identical to explicit setters. */
    ConfigurationPublication publication(*this);
    publication.publish(reason);
}

int
SoBRLDatabaseSource::setDatabaseMetadataState(SbBool metadataValid, int regionId,
    int airCode, int materialId, int los, SbBool colorValid, const SbColor &metadataColor,
    const SbString &shader)
{
    return this->setDatabaseMetadataState(metadataValid, regionId, airCode, materialId,
	los, colorValid, metadataColor, shader, nullptr, nullptr);
}

int
SoBRLDatabaseSource::setDatabaseMetadataState(SbBool metadataValid, int regionId,
    int airCode, int materialId, int los, SbBool colorValid, const SbColor &metadataColor,
    const SbString &shader, PublicationCommit committed, void *context)
{
    const SbColor normalizedColor = colorValid ? metadataColor : SbColor(1.0f, 1.0f, 1.0f);
    const auto visit = [&](auto &node, auto field) {
	field(node.databaseMetadataValid, metadataValid);
	field(node.databaseRegionId, regionId);
	field(node.databaseAirCode, airCode);
	field(node.databaseMaterialId, materialId);
	field(node.databaseLos, los);
	field(node.databaseMaterialColorValid, colorValid);
	field(node.databaseMaterialColor, normalizedColor);
	field(node.databaseMaterialShader, shader);
    };
    bool changed = false;
    visit(*this, [&](const auto &field, const auto &value) {
	changed = publication_value_different(field.getValue(), value) || changed;
    });
    if (!changed) return 0;
    ConfigurationPublication publication(*this);
    visit(publication.next, [](auto &field, const auto &value) {
	if (publication_value_different(field.getValue(), value)) field.setValue(value);
    });
    publication.publishDatabaseMetadata(committed, context);
    return 1;
}

int
SoBRLDatabaseSource::setEvaluatedRegionForPath(const char *queryPath, SbBool evaluated)
{
    return this->setEvaluatedRegionForPath(queryPath, evaluated, nullptr, nullptr);
}

int
SoBRLDatabaseSource::setEvaluatedRegionForPath(const char *queryPath, SbBool evaluated,
    PublicationCommit committed, void *context)
{
    const int region = evaluated ? 1 : 0;
    if (this->hasCompactInstanceIndex())
	return this->setCompactInstanceRegionIdForPath(queryPath, TRUE, region, committed, context);

    struct RegionEdit {
	RegionEdit(SoNode &shape, SoSFInt32 &value) : owner(&shape), field(value),
	    notifications(shape, {{{&value, true}}})
	{}
	SbModernUtils::SoNodeRef owner;
	SoSFInt32 &field;
	PreparedFieldNotifications<1> notifications;
    };
    std::vector<std::unique_ptr<RegionEdit>> changes;
    const auto prepare = [&](auto &shape) {
	if (shape.regionId.getValue() == region)
	    return;
	const SbString &shapePath = shape.sourcePath.getValue().getLength() ?
	    shape.sourcePath.getValue() : this->path.getValue();
	if (database_source_path_matches_frontier(shapePath, queryPath))
	    changes.push_back(std::make_unique<RegionEdit>(shape, shape.regionId));
    };
    std::unordered_set<SoNode *> visited;
    for (int i = 0; i < this->getNumChildren(); ++i) {
	visit_owned_nodes(this->getChild(i), visited, [&](SoNode *owned) {
	    if (node_is_auxiliary_vlist(owned))
		return;
	    if (owned->isOfType(SoBRLVListShape::getClassTypeId()))
		prepare(*static_cast<SoBRLVListShape *>(owned));
	    else if (owned->isOfType(SoBRLMeshShape::getClassTypeId()))
		prepare(*static_cast<SoBRLMeshShape *>(owned));
	});
    }
    if (changes.empty())
	return 0;

    // Only the changed integer fields need prepared storage; scalar metadata
    // copies would duplicate unrelated strings for this one-field operation.
    PreparedFieldNotifications<0> sourceNotifications(*this, {});
    for (auto &change : changes)
	change->field = region;
    this->markCompiledAssemblyDirty();
    this->markCadBatchDirty();
    if (committed)
	committed(context);
    sourceNotifications.restore();
    for (auto &change : changes)
	change->notifications.restore();
    std::exception_ptr failure;
    sourceNotifications.notify(failure);
    for (auto &change : changes)
	change->notifications.notify(failure);
    if (failure)
	std::rethrow_exception(failure);
    return static_cast<int>(changes.size());
}

int
SoBRLDatabaseSource::refreshMaterialColorFromDatabase(uint32_t revision, struct db_i *database)
{
    return this->refreshMaterialColorFromDatabase(revision, database, nullptr, nullptr);
}

int
SoBRLDatabaseSource::refreshMaterialColorFromDatabase(uint32_t revision, struct db_i *database,
    PublicationCommit committed, void *context)
{
    if (this->materialRevision.getValue() == revision) return 0;
    if (!database) database = this->d->dbip;
    if (!database) return 0;
    SoBRLDatabaseSource *source = this;
    return refreshMaterialColors(&source, 1, revision, database, committed, context);
}

int
SoBRLDatabaseSource::refreshMaterialColors(SoBRLDatabaseSource *const *sources, size_t count,
    uint32_t revision, struct db_i *database, PublicationCommit committed, void *context,
    PublicationAccept accepted)
{
    if (!count) return 0;
    if (!sources) return -1;
    bool pending = false;
    for (size_t i = 0; i < count; ++i)
	pending = pending || (sources[i] && sources[i]->materialRevision.getValue() != revision);
    if (!pending) return 0;

    struct Target {
	SbModernUtils::SoNodeRef node;
	uint32_t revision;
	uint64_t cadRevision;
	SbUniqueId nodeId;
    };
    std::vector<Target> targets;
    targets.reserve(count);
    for (size_t i = 0; i < count; ++i) {
	auto *source = sources[i];
	if (source) targets.push_back({SbModernUtils::SoNodeRef(source), source->materialRevision.getValue(),
	    source->cadBatchRevisionGet(), source->getNodeId()});
    }
    struct Progress {
	PublicationCommit committed;
	void *context;
	bool changed = false;
    } progress{committed, context, false};
    const auto published = [](void *data) noexcept {
	auto &state = *static_cast<Progress *>(data);
	state.changed = true;
	if (state.committed) state.committed(state.context);
    };
    std::optional<BObolMaterialColorSweep> sweep;
    struct db_i *sweepDatabase = nullptr;
    for (auto &target : targets) {
	auto &source = *static_cast<SoBRLDatabaseSource *>(target.node.get());
	// A callback can replace a scene owner or make a later quiet edit without
	// changing its material revision. Neither authorizes this older sweep.
	if (source.materialRevision.getValue() == revision || source.materialRevision.getValue() != target.revision ||
	    source.cadBatchRevisionGet() != target.cadRevision || source.getNodeId() != target.nodeId ||
	    (accepted && !accepted(&source, context))) continue;
	auto *targetDatabase = database ? database : source.getDatabase();
	if (!targetDatabase) continue;
	// Reuse one cache while its database remains current. Switching databases
	// retires all path, combination and region-color evidence together.
	if (targetDatabase != sweepDatabase) {
	    sweep.emplace(targetDatabase);
	    sweepDatabase = targetDatabase;
	}
	ConfigurationPublication publication(source);
	publication.refreshMaterials(*sweep, revision, published, &progress);
    }
    return progress.changed ? 1 : 0;
}

int
SoBRLDatabaseSource::setRealizationViewPolicy(SbBool viewDependent, SbBool csgLodEnabled,
    SbBool meshLodEnabled, float viewScale, float lodScale, int viewWidth,
    int viewHeight, uint32_t botThreshold, float curveScale, float pointScale)
{
    return this->setRealizationViewPolicy(viewDependent, csgLodEnabled, meshLodEnabled,
	viewScale, lodScale, viewWidth, viewHeight, botThreshold, curveScale, pointScale,
	nullptr, nullptr);
}

int
SoBRLDatabaseSource::setRealizationViewPolicy(SbBool viewDependent, SbBool csgLodEnabled,
    SbBool meshLodEnabled, float viewScale, float lodScale, int viewWidth,
    int viewHeight, uint32_t botThreshold, float curveScale, float pointScale,
    PublicationCommit committed, void *context)
{
    if (lodScale <= 0.0f) lodScale = 1.0f;
    const auto visit = [&](auto &node, auto field) {
	field(node.realizationViewDependent, viewDependent);
	field(node.realizationCsgLodEnabled, csgLodEnabled);
	field(node.realizationMeshLodEnabled, meshLodEnabled);
	field(node.realizationViewScale, viewScale);
	field(node.realizationLodScale, lodScale);
	field(node.realizationViewWidth, viewWidth);
	field(node.realizationViewHeight, viewHeight);
	field(node.realizationBotThreshold, botThreshold);
	field(node.lodBotThreshold, botThreshold);
	field(node.realizationCurveScale, curveScale);
	field(node.realizationPointScale, pointScale);
    };
    bool changed = false;
    visit(*this, [&](const auto &field, const auto &value) {
	changed = publication_value_different(field.getValue(), value) || changed;
    });
    if (!changed) return 0;
    ConfigurationPublication publication(*this);
    visit(publication.next, [](auto &field, const auto &value) {
	if (publication_value_different(field.getValue(), value)) field.setValue(value);
    });
    publication.publishViewPolicy(committed, context);
    return 1;
}

int
SoBRLDatabaseSource::setDisplayState(SbBool sourceRevisionValid, uint32_t sourceRevisionValue,
    uint32_t inputsRevisionValue, SbBool visibleValue, SbBool selectedValue, SbBool highlightedValue,
    int lineStyleValue, int lineWidthValue, float transparencyValue, SbBool colorOverrideValue,
    const SbColor &colorValue, SbBool materialColorValidValue, const SbColor &materialColorValue,
    uint32_t materialRevisionValue)
{
    return this->setDisplayState(sourceRevisionValid, sourceRevisionValue, inputsRevisionValue,
	visibleValue, selectedValue, highlightedValue, lineStyleValue, lineWidthValue, transparencyValue,
	colorOverrideValue, colorValue, materialColorValidValue, materialColorValue, materialRevisionValue,
	nullptr, nullptr);
}

int
SoBRLDatabaseSource::setDisplayState(SbBool sourceRevisionValid, uint32_t sourceRevisionValue,
    uint32_t inputsRevisionValue, SbBool visibleValue, SbBool selectedValue, SbBool highlightedValue,
    int lineStyleValue, int lineWidthValue, float transparencyValue, SbBool colorOverrideValue,
    const SbColor &colorValue, SbBool materialColorValidValue, const SbColor &materialColorValue,
    uint32_t materialRevisionValue, PublicationCommit committed, void *context)
{
    SourceDisplayChanges changes;
    auto &patch = changes.patch;
    patch.visibleValid = patch.selectedValid = patch.highlightedValid = TRUE;
    patch.lineStyleValid = patch.lineWidthValid = patch.transparencyValid = patch.colorOverrideValid = TRUE;
    patch.visible = visibleValue; patch.selected = selectedValue; patch.highlighted = highlightedValue;
    patch.lineStyle = lineStyleValue; patch.lineWidth = lineWidthValue; patch.transparency = transparencyValue;
    patch.colorOverride = colorOverrideValue; patch.colorValid = colorOverrideValue; patch.color = colorValue;
    const bool displayChanged = filter_display_patch(*this, patch);
    const bool sourceChanged = sourceRevisionValid && this->sourceRevision.getValue() != sourceRevisionValue;
    const bool inputsChanged = this->inputsRevision.getValue() != inputsRevisionValue;
    const bool materialValueChanged = materialColorValidValue &&
	!database_source_color_equal(this->materialColor.getValue(), materialColorValue);
    changes.materialColor = this->materialColorValid.getValue() != materialColorValidValue || materialValueChanged;
    changes.materialRevision = this->materialRevision.getValue() != materialRevisionValue;
    if (!displayChanged && !sourceChanged && !inputsChanged && !changes.materialColor && !changes.materialRevision)
	return 0;
    if (sourceChanged || (this->materialColorValid.getValue() && !materialColorValidValue))
	changes.staleReason |= STALE_SOURCE;
    if (inputsChanged) changes.staleReason |= STALE_INPUTS;

    ConfigurationPublication publication(*this);
    auto &next = publication.next;
    apply_display_patch(next, patch);
    if (sourceChanged) next.sourceRevision = sourceRevisionValue;
    if (inputsChanged) next.inputsRevision = inputsRevisionValue;
    next.materialColorValid = materialColorValidValue;
    if (materialValueChanged) next.materialColor = materialColorValue;
    next.materialRevision = materialRevisionValue;
    return publication.publishDisplay(changes, committed, context) ? 1 : 0;
}

int
SoBRLDatabaseSource::applyDisplayPatch(const BObolDatabaseSourceDisplayPatch &patch)
{
    return this->applyDisplayPatch(patch, nullptr, nullptr);
}

int
SoBRLDatabaseSource::applyDisplayPatch(const BObolDatabaseSourceDisplayPatch &patch,
    PublicationCommit committed, void *context)
{
    SourceDisplayChanges changes;
    changes.patch = patch;
    if (!filter_display_patch(*this, changes.patch)) return 0;
    ConfigurationPublication publication(*this);
    apply_display_patch(publication.next, changes.patch);
    return publication.publishDisplay(changes, committed, context) ? 1 : 0;
}

int
SoBRLDatabaseSource::applyPresentationPatch(
    const BObolDatabaseSourcePresentationPatch &patch,
    PublicationCommit committed, void *context)
{
    if (!presentationPatchValid(patch))
	return -1;

    SourceDisplayChanges changes;
    changes.patch = patch.sourceDisplay;
    (void)filter_display_patch(*this, changes.patch);
    changes.presentation = &patch;

    ConfigurationPublication publication(*this);
    apply_display_patch(publication.next, changes.patch);
    return publication.publishDisplay(changes, committed, context) ? 1 : 0;
}

int
SoBRLDatabaseSource::setDisplayNameState(const char *name)
{
    return this->setDisplayNameState(name, nullptr, nullptr);
}

int
SoBRLDatabaseSource::setDisplayNameState(const char *name, PublicationCommit committed, void *context)
{
    const char *nextName = name ? name : "";
    if (database_source_string_equal(this->displayName.getValue(), nextName))
	return 0;

    ConfigurationPublication publication(*this);
    publication.next.displayName = nextName;
    publication.publishMetadata(MetadataIntent::DisplayName, committed, context);
    return 1;
}

int
SoBRLDatabaseSource::setHierarchyState(
    const char *sourceParentInstanceKey,
    uint32_t sourceOccurrenceIndex,
    int sourceBooleanOperation)
{
    return this->setHierarchyState(sourceParentInstanceKey, sourceOccurrenceIndex, sourceBooleanOperation, nullptr, nullptr);
}

int
SoBRLDatabaseSource::setHierarchyState(const char *sourceParentInstanceKey,
    uint32_t sourceOccurrenceIndex, int sourceBooleanOperation, PublicationCommit committed, void *context)
{
    const char *nextParent = sourceParentInstanceKey ?
	sourceParentInstanceKey : "";
    int nextOperation = sourceBooleanOperation;
    if (nextOperation != BOOLEAN_SUBTRACT &&
	nextOperation != BOOLEAN_INTERSECT)
	nextOperation = BOOLEAN_UNION;

    if (database_source_string_equal(this->parentInstanceKey.getValue(),
	    nextParent) &&
	this->occurrenceIndex.getValue() == sourceOccurrenceIndex &&
	this->booleanOperation.getValue() == nextOperation)
	return 0;

    ConfigurationPublication publication(*this);
    publication.next.parentInstanceKey = nextParent;
    publication.next.occurrenceIndex = sourceOccurrenceIndex;
    publication.next.booleanOperation = nextOperation;
    publication.publishMetadata(MetadataIntent::Hierarchy, committed, context);
    return 1;
}

int
SoBRLDatabaseSource::setMaterialPolicyState(int nextMaterialPolicy)
{
    return this->setMaterialPolicyState(nextMaterialPolicy, nullptr, nullptr);
}

int
SoBRLDatabaseSource::setMaterialPolicyState(int nextMaterialPolicy, PublicationCommit committed, void *context)
{
    if (nextMaterialPolicy != MATERIAL_DATABASE)
	nextMaterialPolicy = MATERIAL_INHERIT;
    if (this->materialPolicy.getValue() == nextMaterialPolicy)
	return 0;

    ConfigurationPublication publication(*this);
    publication.next.materialPolicy = nextMaterialPolicy;
    publication.publishMetadata(MetadataIntent::MaterialPolicy, committed, context);
    return 1;
}

int
SoBRLDatabaseSource::setDrawModeState(int nextDrawMode)
{
    const int mode = nextDrawMode == SHADED ? SHADED : WIREFRAME;
    return this->publishDrawRepresentationState(mode, this->representationKey.getValue().getString(),
	this->representationMode.getValue(), ConfigurationIntent::Draw, nullptr, nullptr);
}

int
SoBRLDatabaseSource::publishDrawRepresentationState(int mode, const char *key, int representation,
    ConfigurationIntent intent, PublicationCommit committed, void *context)
{
    const char *nextKey = key ? key : "";
    const bool drawChanged = this->drawMode.getValue() != mode;
    const bool representationChanged = !database_source_string_equal(this->representationKey.getValue(), nextKey) ||
	this->representationMode.getValue() != representation;
    if (!drawChanged && !representationChanged) return 0;
    const bool preserveExternalRealization = intent == ConfigurationIntent::Draw &&
	source_external_realization_supports_draw_mode(*this, mode);
    ConfigurationPublication publication(*this);
    publication.next.drawMode = mode;
    publication.next.representationKey = nextKey;
    publication.next.representationMode = representation;
    if (preserveExternalRealization)
	publication.next.realizationIdentity = source_realization_identity(&publication.next);
    publication.publish(preserveExternalRealization ? STALE_NONE : STALE_DRAW, representationChanged, committed, context);
    return 1;
}

int
SoBRLDatabaseSource::setRepresentationState(const char *sourceRepresentationKey, int sourceRepresentationMode)
{
    return this->publishDrawRepresentationState(this->drawMode.getValue(), sourceRepresentationKey,
	sourceRepresentationMode, ConfigurationIntent::Representation, nullptr, nullptr);
}

int
SoBRLDatabaseSource::retargetDatabaseSourceInstance(const char *sourceInstanceKey,
    const char *sourcePath, uint32_t revision)
{
    return this->retargetDatabaseSourceInstance(sourceInstanceKey, sourcePath, revision, nullptr, nullptr);
}

int
SoBRLDatabaseSource::retargetDatabaseSourceInstance(const char *sourceInstanceKey,
    const char *sourcePath, uint32_t revision, PublicationCommit committed, void *context)
{
    if (!sourcePath || !sourcePath[0]) return -1;
    const char *instance = sourceInstanceKey && sourceInstanceKey[0] ? sourceInstanceKey : sourcePath;
    const bool pathChanged = !database_source_string_equal(this->path.getValue(), sourcePath);
    const bool instanceChanged = !database_source_string_equal(this->instanceKey.getValue(), instance);
    if (!revision) {
	revision = this->sourceRevision.getValue();
	if (pathChanged || instanceChanged) revision = bobol_identity_successor_or_terminate(revision);
    }
    const bool revisionChanged = this->sourceRevision.getValue() != revision;
    if (!pathChanged && !instanceChanged && !revisionChanged) return 0;
    ConfigurationPublication publication(*this);
    publication.next.instanceKey = instance;
    publication.next.path = sourcePath;
    publication.next.sourceRevision = revision;
    publication.publish(pathChanged || revisionChanged ? STALE_SOURCE : STALE_NONE, true, committed, context);
    return 1;
}

static uint32_t
source_configuration_stale_reason(const SoBRLDatabaseSource &source, const char *instance,
    const char *path, const char *representation, int representationMode, struct db_i *database,
    int drawMode, uint32_t revision)
{
    path = path ? path : "";
    instance = instance && instance[0] ? instance : path;
    representation = representation && representation[0] ? representation : instance;
    drawMode = drawMode == SoBRLDatabaseSource::SHADED ? SoBRLDatabaseSource::SHADED : SoBRLDatabaseSource::WIREFRAME;
    uint32_t reason = SoBRLDatabaseSource::STALE_NONE;
    if (source.getDatabase() != database) reason |= SoBRLDatabaseSource::STALE_DATABASE;
    if (source.drawMode.getValue() != drawMode || source.representationMode.getValue() != representationMode)
	reason |= SoBRLDatabaseSource::STALE_DRAW;
    if (!database_source_string_equal(source.path.getValue(), path) ||
	!database_source_string_equal(source.instanceKey.getValue(), instance) ||
	!database_source_string_equal(source.representationKey.getValue(), representation) ||
	source.sourceRevision.getValue() != revision) reason |= SoBRLDatabaseSource::STALE_SOURCE;
    return reason;
}

void
SoBRLDatabaseSource::configureDatabaseSourceInstanceRepresentation(
    const char *sourceInstanceKey,
    const char *sourcePath,
    const char *sourceRepresentationKey,
    int sourceRepresentationMode,
    struct db_i *database,
    int mode,
    uint32_t revision)
{
    const int sanitizedMode = mode == SHADED ? SHADED : WIREFRAME;
    const char *effectivePath = sourcePath ? sourcePath : "";
    const char *effectiveInstanceKey = sourceInstanceKey && sourceInstanceKey[0] ? sourceInstanceKey : effectivePath;
    const char *effectiveRepresentationKey = sourceRepresentationKey && sourceRepresentationKey[0] ?
	sourceRepresentationKey : effectiveInstanceKey;
    const uint32_t reason = source_configuration_stale_reason(*this, effectiveInstanceKey,
	effectivePath, effectiveRepresentationKey, sourceRepresentationMode, database, sanitizedMode, revision);
    if (!reason && !this->auxiliarySource.getValue()) return;

    ConfigurationPublication publication(*this);
    // Arguments may borrow live fields: copy them all before publishing any value.
    publication.next.instanceKey = effectiveInstanceKey;
    publication.next.path = effectivePath;
    publication.next.representationKey = effectiveRepresentationKey;
    publication.next.representationMode = sourceRepresentationMode;
    publication.next.drawMode = sanitizedMode;
    publication.next.sourceRevision = revision;
    publication.next.auxiliarySource = FALSE;
    publication.database = database;
    publication.publish(reason);
}

/* Placement and child paths are one publication, whether replacing the primary
 * drawing or a named auxiliary. Existing unrelated children remain retained. */
class PreparedSourceChildren {
public:
    PreparedSourceChildren(SoBRLDatabaseSource &source, std::vector<SoNode *> next) :
	PreparedSourceChildren(source, source, std::move(next))
    {}
    PreparedSourceChildren(SoBRLDatabaseSource &source, const SoBRLDatabaseSource &configuration,
	std::vector<SoNode *> next) :
	matrix(configuration.drawMatrix.getValue()), placement(nullptr)
    {
	next.erase(std::remove_if(next.begin(), next.end(), node_is_source_placement_transform), next.end());
	/* Auxiliary source metadata describes the placement already inherited
	 * from its parent; another transform would move its geometry twice. */
	if (configuration.drawMatrixValid.getValue() && !configuration.auxiliarySource.getValue()) {
	    SoMatrixTransform *transform = source_placement_transform(&source);
	    if (!transform) {
		transform = new SoMatrixTransform;
		this->placement = SbModernUtils::SoNodeRef(transform);
		transform->setName(SbName(database_source_placement_transform_name()));
	    } else {
		this->placement = SbModernUtils::SoNodeRef(transform);
	    }
	    this->placementNotifications = std::make_unique<PreparedFieldNotifications<1>>(
		*transform, std::array<PublicationFieldChange, 1>{{{&transform->matrix,
		    (transform->matrix.getValue() != this->matrix) != FALSE}}});
	    next.insert(next.begin(), transform);
	}
	this->nextChildren = std::move(next);
	this->children = source.getChildren()->prepareReplacement(this->nextChildren);
    }
    ~PreparedSourceChildren() { this->restore(); }
    void commit()
    {
	if (this->placement)
	    static_cast<SoMatrixTransform *>(this->placement.get())->matrix = this->matrix;
	this->children->commit();
    }
    void restore()
    {
	if (this->placementNotifications) this->placementNotifications->restore();
    }
    void notify(std::exception_ptr &failure)
    {
	this->restore();
	if (this->placementNotifications) this->placementNotifications->notify(failure);
	try { this->children->notify(); }
	catch (...) { if (!failure) failure = std::current_exception(); }
    }
    void notify()
    {
	std::exception_ptr failure;
	this->notify(failure);
	if (failure) std::rethrow_exception(failure);
    }
    const std::vector<SoNode *> &order() const { return this->nextChildren; }
private:
    SbMatrix matrix;
    SbModernUtils::SoNodeRef placement;
    std::unique_ptr<PreparedFieldNotifications<1>> placementNotifications;
    std::vector<SoNode *> nextChildren;
    std::unique_ptr<SoChildList::Replacement> children;
};

/* Source state, owned metadata and placement commit together. A named auxiliary
 * inherits placement metadata, but never another transform; ordinary nested
 * sources retain their own state. Prepare every participant before committing. */
class SoBRLDatabaseSource::StatePublication {
public:
    explicit StatePublication(SoBRLDatabaseSource &target) : source(target), configuration(target) {}
    SoBRLDatabaseSource &next() { return this->configuration.next; }
    void prepare(std::unordered_set<SoNode *> &visited)
    {
	this->prepare(visited, [&](auto &shape) { sync_shape_placement_state(&shape, &this->next()); }, true);
    }
    template <typename Prepare>
    void prepare(std::unordered_set<SoNode *> &visited, Prepare prepareShape, bool placement)
    {
	std::vector<SoNode *> nextChildren;
	if (placement) nextChildren.reserve(size_t(this->source.getNumChildren()));
	for (int i = 0; i < this->source.getNumChildren(); ++i) {
	    SoNode *child = this->source.getChild(i);
	    if (placement) nextChildren.push_back(child);
	    this->prepareChild(child, visited, prepareShape, placement);
	}
	if (placement)
	    this->children = std::make_unique<PreparedSourceChildren>(this->source, this->next(), std::move(nextChildren));
	this->fields = std::make_unique<PreparedScalarFields>(this->source, this->next());
    }
    void publish(PublicationCommit committed, void *context)
    {
	this->commit();
	if (committed) committed(context);
	this->restore();
	std::exception_ptr failure;
	this->notify(failure);
	if (failure) std::rethrow_exception(failure);
    }
private:
    friend class SoBRLDatabaseSource;
    template <typename Prepare>
    void prepareChild(SoNode *child, std::unordered_set<SoNode *> &visited, Prepare prepareShape, bool placement)
    {
	if (!visited.insert(child).second) return;
	if (child->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    auto &nested = *static_cast<SoBRLDatabaseSource *>(child);
	    if (!placement || !nested.auxiliarySource.getValue()) return;
	    auto edit = std::make_unique<StatePublication>(nested);
	    sync_shape_placement_state(&edit->next(), &this->next());
	    edit->prepare(visited);
	    this->auxiliaries.push_back(std::move(edit));
	} else if (child->isOfType(SoBRLVListShape::getClassTypeId())) {
	    this->prepareShape(*static_cast<SoBRLVListShape *>(child), prepareShape);
	} else if (child->isOfType(SoBRLMeshShape::getClassTypeId())) {
	    this->prepareShape(*static_cast<SoBRLMeshShape *>(child), prepareShape);
	} else if (child->isOfType(SoGroup::getClassTypeId())) {
	    auto &group = *static_cast<SoGroup *>(child);
	    for (int i = 0; i < group.getNumChildren(); ++i) this->prepareChild(group.getChild(i), visited, prepareShape, placement);
	}
    }
    template <typename Shape, typename Prepare>
    void prepareShape(Shape &shape, Prepare prepare)
    {
	auto edit = std::make_unique<PreparedOwnedMetadata>(shape, prepare);
	if (edit->changed()) this->shapes.push_back(std::move(edit));
    }
    void commit()
    {
	this->fields->commit();
	for (auto &shape : this->shapes) shape->commit();
	if (this->children) {
	    this->children->commit();
	    this->source.syncCompactInstancePlacementState();
	    this->source.markDisplayMeshLodDirty();
	}
	this->source.markCompiledAssemblyDirty();
	this->source.markCadBatchDirty();
	for (auto &auxiliary : this->auxiliaries) auxiliary->commit();
    }
    void restore()
    {
	this->fields->restore();
	for (auto &shape : this->shapes) shape->restore();
	if (this->children) this->children->restore();
	for (auto &auxiliary : this->auxiliaries) auxiliary->restore();
    }
    void notify(std::exception_ptr &failure)
    {
	this->fields->notify(failure);
	for (auto &shape : this->shapes) shape->notify(failure);
	if (this->children) this->children->notify(failure);
	for (auto &auxiliary : this->auxiliaries) auxiliary->notify(failure);
    }
    SoBRLDatabaseSource &source;
    ConfigurationPublication configuration;
    std::vector<std::unique_ptr<PreparedOwnedMetadata>> shapes;
    std::vector<std::unique_ptr<StatePublication>> auxiliaries;
    std::unique_ptr<PreparedSourceChildren> children;
    std::unique_ptr<PreparedScalarFields> fields;
};

static bool
source_placement_transform_matches(const SoBRLDatabaseSource &source, SbBool valid, const SbMatrix &matrix)
{
    const bool expected = valid && !source.auxiliarySource.getValue();
    bool found = false;
    for (int i = 0; i < source.getNumChildren(); ++i) {
	const SoNode *child = source.getChild(i);
	if (!node_is_source_placement_transform(child)) continue;
	if (!expected || found || i != 0 ||
	    !static_cast<const SoMatrixTransform *>(child)->matrix.getValue().equals(matrix, 0.000001f))
	    return false;
	found = true;
    }
    return found == expected;
}

static bool
source_placement_state_matches(const SoBRLDatabaseSource &source, SbBool matrixValid,
    const SbMatrix &matrix, SbBool centerValid, const SbVec3f &center, SbBool sizeValid, float size)
{
    return source.drawMatrixValid.getValue() == matrixValid &&
	source.drawMatrix.getValue().equals(matrix, 0.000001f) &&
	source.drawCenterValid.getValue() == centerValid &&
	database_source_vec3f_equal(source.drawCenter.getValue(), center) &&
	source.drawSizeValid.getValue() == sizeValid &&
	!database_source_float_different(source.drawSize.getValue(), size) &&
	source_placement_transform_matches(source, matrixValid, matrix);
}

int
SoBRLDatabaseSource::setPlacementState(SbBool matrixValid, const SbMatrix &matrix,
    SbBool centerValid, const SbVec3f &center, SbBool sizeValid, float size)
{
    return this->setPlacementState(matrixValid, matrix, centerValid, center, sizeValid, size, nullptr, nullptr);
}

int
SoBRLDatabaseSource::setPlacementState(SbBool matrixValid, const SbMatrix &matrix,
    SbBool centerValid, const SbVec3f &center, SbBool sizeValid, float size,
    PublicationCommit committed, void *context)
{
    const bool matrixChanged = !this->drawMatrix.getValue().equals(matrix, 0.000001f);
    const bool centerChanged = !database_source_vec3f_equal(this->drawCenter.getValue(), center);
    const bool sizeChanged = database_source_float_different(this->drawSize.getValue(), size);
    if (source_placement_state_matches(*this, matrixValid, matrix, centerValid, center, sizeValid, size))
	return 0;

    StatePublication publication(*this);
    auto &next = publication.next();
    next.drawMatrixValid = matrixValid;
    if (matrixChanged) next.drawMatrix = matrix;
    next.drawCenterValid = centerValid;
    if (centerChanged) next.drawCenter = center;
    next.drawSizeValid = sizeValid;
    if (sizeChanged) next.drawSize = size;
    std::unordered_set<SoNode *> visited{this};
    publication.prepare(visited);
    publication.publish(committed, context);
    return 1;
}

int
SoBRLDatabaseSource::publishState(const BObolDatabaseSourcePublishState &state, uint32_t revision,
    PublicationCommit committed, void *context)
{
    StatePublication publication(*this);
    auto &next = publication.next();
    (void)next.d->setDatabaseBinding(this->getDatabase());

    // Normalize through the existing setters on a private scalar candidate.
    // No live geometry, compact registry or external observers are copied.
    const bool policyChanged = state.materialPolicyValid && next.setMaterialPolicyState(state.materialPolicy) > 0;
    uint32_t reason = source_configuration_stale_reason(*this, state.sourceInstanceKey,
	state.sourcePath, state.sourceRepresentationKey, state.representationMode,
	state.database, state.drawMode, revision);
    const bool retiresExternalRealization = state.roleFlagsValid &&
	(this->realizationRoleFlags.getValue() & REALIZATION_ROLE_EXTERNAL) &&
	!(state.roleFlags & REALIZATION_ROLE_EXTERNAL);
    if (retiresExternalRealization)
	reason |= STALE_SOURCE;
    const bool configured = reason || this->auxiliarySource.getValue();
    next.configureDatabaseSourceInstanceRepresentation(state.sourceInstanceKey, state.sourcePath,
	state.sourceRepresentationKey, state.representationMode, state.database, state.drawMode, revision);
    if (retiresExternalRealization)
	next.markStale(STALE_SOURCE);

    bool displayChanged = false;
    if (state.sourceRevisionValid || state.inputsRevision != 0 ||
	next.visible.getValue() != state.visible || next.selected.getValue() != state.selected ||
	next.highlighted.getValue() != state.highlighted || next.lineStyle.getValue() != state.lineStyle ||
	next.lineWidth.getValue() != state.lineWidth ||
	database_source_float_different(next.transparency.getValue(), state.transparency) ||
	next.colorOverride.getValue() != state.colorOverride ||
	(state.colorOverride && !database_source_color_equal(next.color.getValue(), state.color)) ||
	next.materialColorValid.getValue() != state.materialColorValid ||
	(state.materialColorValid && !database_source_color_equal(next.materialColor.getValue(), state.materialColor)) ||
	next.materialRevision.getValue() != state.materialRevision) {
	if (next.inputsRevision.getValue() != state.inputsRevision) reason |= STALE_INPUTS;
	if (next.materialColorValid.getValue() && !state.materialColorValid) reason |= STALE_SOURCE;
	displayChanged = next.setDisplayState(state.sourceRevisionValid, revision, state.inputsRevision,
	    state.visible, state.selected, state.highlighted, state.lineStyle, state.lineWidth,
	    state.transparency, state.colorOverride, state.color, state.materialColorValid,
	    state.materialColor, state.materialRevision) > 0;
    }
    const bool rolesChanged = state.roleFlagsValid && next.setRealizationRoleFlags(state.roleFlags) > 0;
    const bool viewChanged = state.viewPolicyValid && next.setRealizationViewPolicy(state.viewDependent,
	state.csgLodEnabled, state.meshLodEnabled, state.viewScale, state.lodScale, state.viewWidth,
	state.viewHeight, state.botThreshold, state.curveScale, state.pointScale) > 0;
    const bool realizationChanged = state.realizationStateValid &&
	next.setRealizationState(state.realizationStatus,
	    state.realizedSourceRevision, state.realizedInputsRevision,
	    state.staleReason, state.realizationDiagnostic) > 0;
    const bool placementChanged = state.placementValid &&
	!source_placement_state_matches(*this, state.drawMatrixValid, state.drawMatrix,
	    state.drawCenterValid, state.drawCenter, state.drawSizeValid, state.drawSize);
    if (placementChanged)
	next.setPlacementState(state.drawMatrixValid, state.drawMatrix, state.drawCenterValid,
	    state.drawCenter, state.drawSizeValid, state.drawSize);
    if (!configured && !policyChanged && !displayChanged && !rolesChanged &&
	!viewChanged && !realizationChanged && !placementChanged)
	return 0;

    SourceDisplayChanges display;
    visit_display_patch(next, display.patch, [](const auto &field, SbBool &valid, auto &value) {
	valid = TRUE;
	value = field.getValue();
    });
    (void)filter_display_patch(*this, display.patch);
    display.materialColor = this->materialColorValid.getValue() != next.materialColorValid.getValue() ||
	!database_source_color_equal(this->materialColor.getValue(), next.materialColor.getValue());
    display.materialRevision = this->materialRevision.getValue() != next.materialRevision.getValue();
    BObolMaterialColorSweep sweep(state.database);
    std::unordered_set<SoNode *> visited{this};
    publication.prepare(visited, [&](auto &shape) {
	if (configured) sync_shape_source_state(&shape, &next);
	else if (reason || viewChanged || realizationChanged)
	    sync_shape_owner_state(&shape, &next);
	if (placementChanged) sync_shape_placement_state(&shape, &next);
	if (!shape.databaseIntent.getValue() && shape.nonDatabaseSource.getValue()) return;
	apply_display_patch(shape, display.patch);
	if (policyChanged || display.materialColor) prepare_shape_material_policy(shape, next, sweep);
	else if (display.materialRevision) shape.materialRevision = next.materialRevision.getValue();
    }, placementChanged);
    PreparedCompactDisplay compact(this->d->compactIndex, next, display, *this->d, this->isCompactOccurrenceRegistry());

    (void)this->d->setDatabaseBinding(state.database);
    publication.commit();
    if (configured) compact_sync_draw_mode(this->d->compactIndex, *this);
    compact.commit();
    if (policyChanged) compact_sync_material_policy(this->d->compactIndex, *this);
    if (invalidation_revokes_source_contract(reason))
	this->revokeSourceResourceAuthority();
    this->markDisplayMeshLodVisibilityDirty(compact.visibilityChangedEntries());
    if (committed) committed(context);
    publication.restore();
    std::exception_ptr failure;
    publication.notify(failure);
    if (failure) std::rethrow_exception(failure);
    return 1;
}

class BObolPreparedSourcePublication {
public:
    BObolPreparedSourcePublication(SoBRLDatabaseSource &target,
	const std::vector<SoNode *> &primary, const SbBox3f &bounds,
	bool preserveAuxiliary = true) :
	BObolPreparedSourcePublication(target, target, primary, bounds, TRUE,
	    SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL, false, preserveAuxiliary)
    {}

private:
    BObolPreparedSourcePublication(SoBRLDatabaseSource &target,
	const SoBRLDatabaseSource &next, const std::vector<SoNode *> &primary, const SbBox3f &bounds,
	SbBool exact, int roles, bool clearSelection, bool preserveAuxiliary) :
	source(target), realization(target, next.sourceRevision.getValue(), bounds, exact, clearSelection, roles)
    {
	std::vector<SoNode *> nextChildren;
	std::unordered_set<SoNode *> preparedAuxiliary;
	for (int i = 0; preserveAuxiliary && i < target.getNumChildren(); ++i) {
	    SoNode *child = target.getChild(i);
	    if (node_is_auxiliary_source(child)) {
		nextChildren.push_back(child);
	    } else if (node_is_auxiliary_vlist(child)) {
		nextChildren.push_back(child);
		if (preparedAuxiliary.insert(child).second) {
		    auto metadata = std::make_unique<PreparedOwnedMetadata>(
			*static_cast<SoBRLVListShape *>(child), target, next.sourceRevision.getValue(), false);
		    if (metadata->changed()) this->changedNodes.push_back(child);
		    this->auxiliary.push_back(std::move(metadata));
		}
	    }
	}
	nextChildren.insert(nextChildren.end(), primary.begin(), primary.end());
	this->children = std::make_unique<PreparedSourceChildren>(target, std::move(nextChildren));
    }
public:
    ~BObolPreparedSourcePublication() { this->restore(); }

    template <typename Construct>
    static int realize(SoBRLDatabaseSource *target, Construct construct,
	BObolSourceRealizationEffects *effects = nullptr)
    {
	if (!target)
	    return -1;
	auto owner = prepareCandidate(*target);
	auto *detached = static_cast<SoBRLDatabaseSource *>(owner.get());
	const int result = construct(detached);
	if (result <= 0) {
	    if (effects) {
		const std::vector<SoNode *> none;
		effects->prepare(*target, false, detached->realizationDiagnostic.getValue(), none, none);
	    }
	    const int changed = target->publishRealizationState(SoBRLDatabaseSource::FAILED, 0, 0,
		target->staleReason.getValue(), detached->realizationDiagnostic.getValue().getString(), effects);
	    if (effects && !changed) effects->commit(false);
	    return result;
	}
	publish(*target, *detached, effects);
	return result;
    }

    static int refreshCombination(SoBRLDatabaseSource &target, uint32_t revision)
    {
	BObolDatabaseSourceRealizationCache cache;
	target.seedCompactRealizationCache(&cache);
	auto owner = prepareCandidate(target);
	auto &detached = *static_cast<SoBRLDatabaseSource *>(owner.get());
	detached.sourceRevision = revision;
	/* The old exact bound describes the preceding combination placement. */
	detached.clearSourceBounds();
	const int result = target.usesMeshRealization() ?
	    bobol_database_source_construct_mesh_compact_with_cache(&detached, &cache, nullptr) :
	    bobol_database_source_construct_wireframe_compact_with_cache(&detached, &cache, nullptr);
	/* A failed synchronous edit keeps its preceding drawing and revision.
	 * The caller owns subsequent invalidation or fallback realization. */
	if (result <= 0)
	    return -1;
	publish(target, detached);
	return target.getCompactInstanceCount();
    }

    static int publishCompactSnapshot(SoBRLDatabaseSource &target,
	const BObolSourceRealizationStamp &stamp,
	const std::vector<BObolCompactOccurrence> &occurrences,
	const SbBox3f *certifiedBounds,
	const BObolCompactSourceProfile *profile,
	BObolSourceRealizationEffects *effects)
    {
	if (occurrences.empty() || !target.matchesRealizationStamp(stamp))
	    return 0;
	auto owner = prepareCandidate(target);
	auto &detached = *static_cast<SoBRLDatabaseSource *>(owner.get());
	const int result = detached.setCompactOccurrenceRegistry(occurrences);
	if (result <= 0)
	    return result;
	if (certifiedBounds) {
	    (void)detached.setSourceBoundsState(TRUE,
		certifiedBounds->getMin(), certifiedBounds->getMax(), TRUE);
	}
	if (profile && detached.certifyCompactStream(
		static_cast<size_t>(profile->occurrenceCount), profile) < 0)
	    return 0;
	detached.realizationRoleFlags =
	    target.realizationRoleFlags.getValue() |
	    SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL;
	/* Preparation is private and normally cannot dispatch a scene edit. Keep
	 * the admission check at the publication boundary so future constructors
	 * cannot accidentally weaken stamped acceptance. */
	if (!target.matchesRealizationStamp(stamp))
	    return 0;
	publish(target, detached, effects);
	return result;
    }

private:
    static SbModernUtils::SoNodeRef prepareCandidate(SoBRLDatabaseSource &target)
    {
	/* Existing geometry stays live; only configuration is copied. This
	 * owning-thread operation borrows the database for private construction. */
	SoBRLDatabaseSource *detached = target.createDetachedRealizationTemplate();
	SbModernUtils::SoNodeRef owner(detached);
	detached->unref();
	(void)detached->d->setDatabaseBinding(target.getDatabase());
	return owner;
    }

    static void publish(SoBRLDatabaseSource &target, SoBRLDatabaseSource &detached,
	BObolSourceRealizationEffects *effects = nullptr)
    {
	std::vector<SoNode *> primary;
	for (int i = 0; i < detached.getNumChildren(); ++i) {
	    SoNode *child = detached.getChild(i);
	    if (!node_is_source_placement_transform(child))
		primary.push_back(child);
	}
	SbBox3f bounds;
	(void)detached.getSourceBounds(bounds);
	BObolPreparedSourcePublication publication(target, detached, primary, bounds,
	    detached.hasExactSourceBounds(), detached.realizationRoleFlags.getValue(),
	    detached.d->compactIndex && !target.d->compactSelectedPaths.empty(), true);
	if (effects) effects->prepare(target, true, detached.realizationDiagnostic.getValue(),
	    publication.children->order(), publication.changedNodes);
	std::unique_ptr<BObolCompactInstanceIndex> index(detached.d->compactIndex);
	detached.d->compactIndex = nullptr;
	if (index) {
	    /* Installation prepares overrides and stable handles against the live
	     * registry before its first write. Remaining publication is quiet and
	     * does not allocate, including child/path and compiled retirement. */
	    target.installCompactInstanceIndex(std::move(index),
		detached.d->compactOccurrenceRegistry, detached.d->compactHandleSourceId);
	} else {
	    target.discardCompactInstanceHistory();
	}
	/* Installation binds the source-mesh contract to the published revision,
	 * which an edit may advance along with its new geometry. */
	target.d->displayMeshLodContractSourceRevision = detached.sourceRevision.getValue();
	target.d->compactExpectedInstanceCount = detached.d->compactExpectedInstanceCount;
	target.d->compactExpectedInstanceCountCertified = detached.d->compactExpectedInstanceCountCertified;
	target.d->compactSourceProfile = detached.d->compactSourceProfile;
	publication.commitFieldsAndChildren();
	if (effects) effects->commit(true);
	std::exception_ptr failure;
	publication.notify(failure);
	if (effects) effects->notify(failure);
	if (failure) std::rethrow_exception(failure);
    }

public:
    void stageEffects(BObolSourceChildEffects *effects) const
    {
	if (!effects) return;
	effects->stageChildOrder(this->source, this->children->order());
	for (auto *node : this->changedNodes) effects->stageFrameEffect(*node);
    }
    void commit()
    {
	this->source.discardCompactInstanceHistory();
	this->commitFieldsAndChildren();
    }
private:
    void commitFieldsAndChildren()
    {
	this->realization.commit();
	for (auto &shape : this->auxiliary) shape->commit();
	this->children->commit();
	/* Retire compiled children only after the prepared child edit removed them,
	 * so cleanup cannot dispatch a legacy partial-removal notification. */
	this->source.clearCompiledAssembly();
    }
public:
    void restore()
    {
	this->realization.restore();
	if (this->children) this->children->restore();
	for (auto &shape : this->auxiliary) shape->restore();
    }
    void notify(std::exception_ptr &failure)
    {
	this->restore();
	this->realization.notify(failure);
	for (auto &shape : this->auxiliary) shape->notify(failure);
	this->children->notify(failure);
    }
    void notify()
    {
	std::exception_ptr failure;
	this->notify(failure);
	if (failure) std::rethrow_exception(failure);
    }
private:
    SoBRLDatabaseSource &source;
    PreparedSourceRealization realization;
    std::vector<std::unique_ptr<PreparedOwnedMetadata>> auxiliary;
    std::vector<SoNode *> changedNodes;
    std::unique_ptr<PreparedSourceChildren> children;
};

int
SoBRLDatabaseSource::publishCompactSnapshot(
    const BObolSourceRealizationStamp &stamp,
    const std::vector<BObolCompactOccurrence> &occurrences,
    const SbBox3f *certifiedBounds,
    const BObolCompactSourceProfile *profile,
    BObolSourceRealizationEffects *effects)
{
    return BObolPreparedSourcePublication::publishCompactSnapshot(*this,
	stamp, occurrences, certifiedBounds, profile, effects);
}

SbBool
SoBRLDatabaseSource::realizePrototypeWireframe(void)
{
    return bobol_database_source_realize_prototype(this);
}

SbBool
bobol_database_source_realize_prototype(SoBRLDatabaseSource *source,
    BObolSourceRealizationEffects *effects)
{
    return BObolPreparedSourcePublication::realize(source,
	construct_prototype_wireframe, effects) > 0 ? TRUE : FALSE;
}

SbBool
bobol_database_source_realize_wireframe_with_cache(
    SoBRLDatabaseSource *source, BObolDatabaseSourceRealizationCache *cache,
    BObolSourceRealizationEffects *effects)
{
    return BObolPreparedSourcePublication::realize(source,
	[&](SoBRLDatabaseSource *detached) {
	    return bobol_database_source_construct_wireframe_with_cache(detached, cache);
	}, effects) > 0 ? TRUE : FALSE;
}

SbBool
bobol_database_source_realize_mesh_with_cache(
    SoBRLDatabaseSource *source, BObolDatabaseSourceRealizationCache *cache,
    BObolSourceRealizationEffects *effects)
{
    return BObolPreparedSourcePublication::realize(source,
	[&](SoBRLDatabaseSource *detached) {
	    return bobol_database_source_construct_mesh_with_cache(detached, cache);
	}, effects) > 0 ? TRUE : FALSE;
}

int
bobol_database_source_realize_wireframe_compact_with_cache(
    SoBRLDatabaseSource *source, BObolDatabaseSourceRealizationCache *cache,
    BObolCompactOccurrenceStream *stream, BObolSourceRealizationEffects *effects)
{
    if (source && source_uses_evaluated_path_realization(source))
	return 0;
    return BObolPreparedSourcePublication::realize(source,
	[&](SoBRLDatabaseSource *detached) {
	    return bobol_database_source_construct_wireframe_compact_with_cache(detached, cache, stream);
	}, effects);
}

int
bobol_database_source_realize_mesh_compact_with_cache(
    SoBRLDatabaseSource *source, BObolDatabaseSourceRealizationCache *cache,
    BObolCompactOccurrenceStream *stream, BObolSourceRealizationEffects *effects)
{
    if (source && source_uses_evaluated_path_realization(source))
	return 0;
    return BObolPreparedSourcePublication::realize(source,
	[&](SoBRLDatabaseSource *detached) {
	    return bobol_database_source_construct_mesh_compact_with_cache(detached, cache, stream);
	}, effects);
}

SbBool
bobol_database_source_construct_realization(SoBRLDatabaseSource *source,
    SbBool mesh, BObolCompactOccurrenceStream *stream)
{
    BObolDatabaseSourceRealizationCache cache;
    if (mesh) {
	if (source_uses_evaluated_points_realization(source))
	    return bobol_database_source_construct_mesh_with_cache(source, &cache);
	return bobol_database_source_construct_mesh_compact_with_cache(source, &cache, stream) > 0;
    }
    if (source_uses_evaluated_wire_realization(source))
	return bobol_database_source_construct_wireframe_with_cache(source, &cache);
    return bobol_database_source_construct_wireframe_compact_with_cache(source, &cache, stream) > 0;
}

static SbBox3f
external_bounds(const SbVec3f *points, int count)
{
    SbVec3f minimum, maximum;
    return external_bounds_from_points(points, count, minimum, maximum) ?
	SbBox3f(minimum, maximum) : SbBox3f();
}

template <typename Shape>
static int
publish_prepared_external_shape(SoBRLDatabaseSource &source, Shape &shape,
    const SbBox3f &bounds, BObolSourceChildEffects *effects)
{
    prepare_realized_shape_owner(shape, source);
    BObolPreparedSourcePublication publication(source, {&shape}, bounds);
    publication.stageEffects(effects);
    if (effects) effects->prepare();
    publication.commit();
    if (effects) effects->commit();
    publication.notify();
    return 1;
}

static SoBRLVListShape *
first_direct_primary_vlist_child(SoBRLDatabaseSource *source)
{
    if (!source)
	return NULL;

    for (int i = 0; i < source->getNumChildren(); i++) {
	SoNode *child = source->getChild(i);
	if (!child || !child->isOfType(SoBRLVListShape::getClassTypeId()))
	    continue;
	if (node_is_auxiliary_vlist(child))
	    continue;
	return static_cast<SoBRLVListShape *>(child);
    }

    return NULL;
}

static SoBRLMeshShape *
first_direct_primary_mesh_child(SoBRLDatabaseSource *source)
{
    if (!source)
	return NULL;

    for (int i = 0; i < source->getNumChildren(); i++) {
	SoNode *child = source->getChild(i);
	if (child && child->isOfType(SoBRLMeshShape::getClassTypeId()))
	    return static_cast<SoBRLMeshShape *>(child);
    }

    return NULL;
}

BObolSceneLightRealization::BObolSceneLightRealization(void) :
    kind(BOBOL_SCENE_LIGHT_POINT),
    position(0.0f, 0.0f, 0.0f),
    direction(0.0f, 0.0f, -1.0f),
    color(1.0f, 1.0f, 1.0f),
    intensity(1.0f),
    coneAngleDeg(180.0f)
{
}

int
SoBRLDatabaseSource::clearRealizedGeometry(SbBool preserveAuxiliary)
{
    std::vector<int> removedIndices;
    int removedRealizedGeometry = 0;
    for (int i = 0; i < this->getNumChildren(); ++i) {
	SoNode *child = this->getChild(i);
	const bool retained = node_is_source_placement_transform(child) ||
	    (preserveAuxiliary && source_retains_realized_child(child));
	if (retained)
	    continue;
	removedIndices.push_back(i);
	if (child != this->d->compiledAssembly)
	    removedRealizedGeometry = 1;
    }

    /* Prepare every path and lifetime change before touching live state.  The
     * compact registry and compiled assembly do not allocate while clearing,
     * so observers can run only after the complete retirement is visible. */
    auto removal = this->getChildren()->prepareRemoval(std::move(removedIndices));
    if (removal)
	removal->commit();
    this->discardCompactInstanceHistory();
    this->clearCompiledAssembly();
    if (removal)
	removal->notify();
    return removedRealizedGeometry;
}

class PreparedExternalPrimaryClear {
public:
    explicit PreparedExternalPrimaryClear(SoBRLDatabaseSource &source) : empty(nullptr)
    {
	SoBRLVListShape *wire = first_direct_primary_vlist_child(&source);
	SoBRLMeshShape *mesh = first_direct_primary_mesh_child(&source);
	if (!wire && !mesh) {
	    wire = new SoBRLVListShape;
	    this->empty = SbModernUtils::SoNodeRef(wire);
	    assign_external_primary_identity(wire, &source, "line-set", "line");
	    prepare_realized_shape_owner(*wire, source);
	    this->nodes.push_back(wire);
	    return;
	}
	if (wire) {
	    this->nodes.push_back(wire);
	    this->clearing.push_back(std::make_unique<PreparedOwnedMetadata>(*wire, source, true));
	}
	if (mesh) {
	    this->nodes.push_back(mesh);
	    this->clearing.push_back(std::make_unique<PreparedOwnedMetadata>(*mesh, source, true));
	}
    }
    void publish(BObolPreparedSourcePublication &publication,
	BObolSourceChildEffects *effects)
    {
	if (effects) {
	    for (const auto &shape : this->clearing)
		if (shape->changed()) effects->stageFrameEffect(*shape->target());
	    publication.stageEffects(effects);
	    effects->prepare();
	}
	for (auto &shape : this->clearing) shape->commit();
	publication.commit();
	if (effects) effects->commit();
	for (auto &shape : this->clearing) shape->restore();
	std::exception_ptr failure;
	publication.notify(failure);
	for (auto &shape : this->clearing) shape->notify(failure);
	if (failure) std::rethrow_exception(failure);
    }
    std::vector<SoNode *> nodes;
private:
    SbModernUtils::SoNodeRef empty;
    std::vector<std::unique_ptr<PreparedOwnedMetadata>> clearing;
};

static int
publish_external_primary_clear(SoBRLDatabaseSource &source,
    BObolSourceChildEffects *effects)
{
    PreparedExternalPrimaryClear empty(source);
    BObolPreparedSourcePublication publication(source, empty.nodes, SbBox3f());
    empty.publish(publication, effects);
    return 1;
}

int
SoBRLDatabaseSource::clearExternalPrimaryGeometry(void)
{
    return this->clearExternalPrimaryGeometry(nullptr);
}

int
SoBRLDatabaseSource::clearExternalPrimaryGeometry(
    BObolSourceChildEffects *effects)
{
    return publish_external_primary_clear(*this, effects);
}

int
SoBRLDatabaseSource::publishExternalLineSet(
    const BObolExternalLineSet &lineSet)
{
    return this->publishExternalLineSet(lineSet, nullptr);
}

int
SoBRLDatabaseSource::publishExternalLineSet(
    const BObolExternalLineSet &lineSet,
    BObolSourceChildEffects *effects)
{
    if (lineSet.count < 0 || (lineSet.count > 0 && !lineSet.points))
	return 0;

    std::vector<int32_t> fallbackCommands;
    const int32_t *commands = lineSet.commands;
    if (lineSet.count > 0) {
	if (!commands) {
	    fallbackCommands.reserve(lineSet.count);
	    for (int i = 0; i < lineSet.count; i++)
		fallbackCommands.push_back(i == 0 ? SoBRLVListShape::MOVE :
					   SoBRLVListShape::DRAW);
	    commands = fallbackCommands.data();
	}
	for (int i = 0; i < lineSet.count; i++) {
	    if (!external_vlist_command_valid(commands[i]))
		return 0;
	}
    }

    if (lineSet.count == 0)
	return this->clearExternalPrimaryGeometry(effects);

    SoBRLVListShape *shape = new SoBRLVListShape;
    SbModernUtils::SoNodeRef preparedShape(shape);

    assign_external_primary_identity(shape, this,
				     external_string_or_default(lineSet.sourceType, "line-set"),
				     external_string_or_default(lineSet.geometryKind, "line"));
    shape->setLineSet(lineSet.points, commands, lineSet.count);
    shape->setPrecisePoints(lineSet.precisePoints, lineSet.count);
    return publish_prepared_external_shape(*this, *shape,
	external_bounds(lineSet.points, lineSet.count), effects);
}

int
SoBRLDatabaseSource::publishExternalPointSet(
    const BObolExternalPointSet &pointSet)
{
    return this->publishExternalPointSet(pointSet, nullptr);
}

int
SoBRLDatabaseSource::publishExternalPointSet(
    const BObolExternalPointSet &pointSet,
    BObolSourceChildEffects *effects)
{
    if (pointSet.count < 0 || (pointSet.count > 0 && !pointSet.points))
	return 0;

    if (pointSet.count == 0)
	return this->clearExternalPrimaryGeometry(effects);

    std::vector<int32_t> commands;
    commands.reserve(pointSet.count);
    for (int i = 0; i < pointSet.count; i++)
	commands.push_back(SoBRLVListShape::POINT);

    BObolExternalLineSet lineSet;
    lineSet.points = pointSet.points;
    lineSet.commands = commands.data();
    lineSet.precisePoints = pointSet.precisePoints;
    lineSet.count = pointSet.count;
    lineSet.sourceType =
	external_string_or_default(pointSet.sourceType, "point-set");
    lineSet.geometryKind =
	external_string_or_default(pointSet.geometryKind, "point");
    return this->publishExternalLineSet(lineSet, effects);
}

int
SoBRLDatabaseSource::publishExternalTriangleMesh(
    const BObolExternalTriangleMesh &triangleMesh)
{
    return this->publishExternalTriangleMesh(triangleMesh, nullptr);
}

int
SoBRLDatabaseSource::publishExternalTriangleMesh(
    const BObolExternalTriangleMesh &triangleMesh,
    BObolSourceChildEffects *effects)
{
    if (triangleMesh.pointCount < 0 || triangleMesh.indexCount < 0 ||
	triangleMesh.normalCount < 0 ||
	(triangleMesh.pointCount > 0 && !triangleMesh.points) ||
	(triangleMesh.indexCount > 0 && !triangleMesh.indices) ||
	(triangleMesh.normalCount > 0 && !triangleMesh.normals) ||
	(triangleMesh.normalCount > 0 &&
	 triangleMesh.normalCount != triangleMesh.indexCount) ||
	(triangleMesh.indexCount % 3) != 0)
	return 0;

    if (triangleMesh.pointCount == 0 || triangleMesh.indexCount == 0)
	return this->clearExternalPrimaryGeometry(effects);

    for (int i = 0; i < triangleMesh.indexCount; i++) {
	const int32_t index = triangleMesh.indices[i];
	if (index < 0 || index >= triangleMesh.pointCount)
	    return 0;
    }

    SoBRLMeshShape *shape = triangleMesh.lodBacked ?
			    new SoBRLLodMeshShape : new SoBRLMeshShape;
    SbModernUtils::SoNodeRef preparedShape(shape);
    if (triangleMesh.lodBacked)
	shape->setLodBackedMesh(TRUE);

    assign_external_primary_identity(shape, this,
				     external_string_or_default(triangleMesh.sourceType,
					     "indexed-face-set"),
				     external_string_or_default(triangleMesh.geometryKind, "surface"));
    std::vector<SbVec3f> points(triangleMesh.points,
	triangleMesh.points + triangleMesh.pointCount);
    std::vector<int32_t> indices(triangleMesh.indices,
	triangleMesh.indices + triangleMesh.indexCount);
    std::vector<SbVec3f> normals;
    if (triangleMesh.normalCount > 0)
	normals.assign(triangleMesh.normals,
	    triangleMesh.normals + triangleMesh.normalCount);
    sanitize_triangle_normals(normals, points, indices);
    shape->setIndexedTriangles(points.data(), static_cast<int>(points.size()),
	indices.data(), static_cast<int>(indices.size()),
	normals.empty() ? NULL : normals.data(),
	static_cast<int>(normals.size()));
    return publish_prepared_external_shape(*this, *shape,
	external_bounds(triangleMesh.points, triangleMesh.pointCount), effects);
}

int
SoBRLDatabaseSource::publishExternalAnnotation(
    const BObolExternalAnnotation &annotation)
{
    return this->publishExternalAnnotation(annotation, nullptr);
}

int
SoBRLDatabaseSource::publishExternalAnnotation(
    const BObolExternalAnnotation &annotation,
    BObolSourceChildEffects *effects)
{
    if (annotation.linePointCount < 0 || annotation.annotationPointCount < 0 ||
	annotation.segmentCount < 0 ||
	(annotation.linePointCount > 0 && !annotation.linePoints) ||
	(annotation.annotationPointCount > 0 &&
	 !annotation.annotationPoints) ||
	(annotation.segmentCount > 0 && !annotation.segments))
	return 0;

    if (annotation.lineCommands) {
	for (int i = 0; i < annotation.linePointCount; i++) {
	    if (!external_vlist_command_valid(annotation.lineCommands[i]))
		return 0;
	}
    }

    std::vector<int32_t> fallbackCommands;
    const int32_t *lineCommands = annotation.lineCommands;
    if (annotation.linePointCount > 0 && !lineCommands) {
	fallbackCommands.reserve(annotation.linePointCount);
	for (int i = 0; i < annotation.linePointCount; i++)
	    fallbackCommands.push_back(i == 0 ? SoBRLVListShape::MOVE :
				       SoBRLVListShape::DRAW);
	lineCommands = fallbackCommands.data();
    }

    SoBRLVListShape *shape = new SoBRLVListShape;
    SbModernUtils::SoNodeRef preparedShape(shape);

    assign_external_primary_identity(shape, this,
				     external_string_or_default(annotation.sourceType, "annotation"),
				     external_string_or_default(annotation.geometryKind, "annotation"));
    shape->setLineSet(annotation.linePoints, lineCommands,
		      annotation.linePointCount);
    shape->setPrecisePoints(annotation.preciseLinePoints,
			    annotation.linePointCount);
    shape->annotationBasePoint = annotation.basePoint;
    if (annotation.annotationPointCount > 0)
	shape->annotationPoint.setValues(0, annotation.annotationPointCount,
					 annotation.annotationPoints);
    else
	shape->annotationPoint.setNum(0);
    if (annotation.preciseAnnotationPoints)
	shape->setPreciseAnnotationPoints(annotation.preciseAnnotationPoints,
					  annotation.annotationPointCount);

    shape->annotationSegmentKind.setNum(annotation.segmentCount);
    shape->annotationSegmentStart.setNum(annotation.segmentCount);
    shape->annotationSegmentEnd.setNum(annotation.segmentCount);
    shape->annotationTextRefPoint.setNum(annotation.segmentCount);
    shape->annotationText.setNum(annotation.segmentCount);
    shape->annotationSegmentTextValid.setNum(annotation.segmentCount);
    for (int i = 0; i < annotation.segmentCount; i++) {
	const BObolExternalAnnotationSegment &segment =
	    annotation.segments[i];
	int kind = SoBRLVListShape::ANNOTATION_SEGMENT_NONE;
	if (segment.kind == BObolExternalAnnotationSegment::SEGMENT_LINE)
	    kind = SoBRLVListShape::ANNOTATION_SEGMENT_LINE;
	else if (segment.kind == BObolExternalAnnotationSegment::SEGMENT_TEXT)
	    kind = SoBRLVListShape::ANNOTATION_SEGMENT_TEXT;
	shape->annotationSegmentKind.set1Value(i, kind);
	shape->annotationSegmentStart.set1Value(i, segment.lineStart);
	shape->annotationSegmentEnd.set1Value(i, segment.lineEnd);
	shape->annotationTextRefPoint.set1Value(i, segment.textRefPoint);
	shape->annotationText.set1Value(i,
					(segment.text && segment.text[0]) ? segment.text : "");
	shape->annotationSegmentTextValid.set1Value(i,
		(kind == SoBRLVListShape::ANNOTATION_SEGMENT_TEXT &&
		 segment.text && segment.text[0]) ? TRUE : FALSE);
    }

    return publish_prepared_external_shape(*this, *shape,
	external_bounds(annotation.linePoints, annotation.linePointCount), effects);
}

static void
primitive_realization_line_set_free(
    struct rt_primitive_lod_realization *realization)
{
    rt_primitive_lod_realization_free(realization);
}

static int32_t
primitive_realization_command_to_vlist_command(int command)
{
    switch (command) {
	case RT_PRIMITIVE_LINE_MOVE:
	    return SoBRLVListShape::MOVE;
	case RT_PRIMITIVE_LINE_DRAW:
	    return SoBRLVListShape::DRAW;
	case RT_PRIMITIVE_POINT_DRAW:
	    return SoBRLVListShape::POINT;
	default:
	    break;
    }
    return -1;
}

static SoBRLVListShape *
vlist_from_primitive_realization_line_set(
    struct rt_primitive_lod_realization *realization,
    const char *geometryKind)
{
    if (!realization || !realization->has_line_set)
	return NULL;

    std::unique_ptr<rt_primitive_lod_realization, decltype(&primitive_realization_line_set_free)> storage(
	realization, primitive_realization_line_set_free);
    if (realization->line_count > static_cast<size_t>(INT_MAX)) {
	return NULL;
    }

    if (realization->line_count == 0 || !realization->line_points) {
	return NULL;
    }

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    std::vector<double> precisePoints;
    points.reserve(realization->line_count);
    commands.reserve(realization->line_count);
    precisePoints.reserve(realization->line_count * 3);
    for (size_t i = 0; i < realization->line_count; i++) {
	const int32_t command =
	    primitive_realization_command_to_vlist_command(
		realization->line_commands ? realization->line_commands[i] :
		RT_PRIMITIVE_LINE_DRAW);
	if (command < 0) {
	    return NULL;
	}
	points.push_back(SbVec3f(
			     static_cast<float>(realization->line_points[i][X]),
			     static_cast<float>(realization->line_points[i][Y]),
			     static_cast<float>(realization->line_points[i][Z])));
	commands.push_back(command);
	precisePoints.push_back(realization->line_points[i][X]);
	precisePoints.push_back(realization->line_points[i][Y]);
	precisePoints.push_back(realization->line_points[i][Z]);
    }

    SoBRLVListShape *shape = new SoBRLVListShape;
    SbModernUtils::SoNodeRef prepared(shape);
    shape->setLineSet(points.empty() ? NULL : points.data(),
		      commands.empty() ? NULL : commands.data(),
		      static_cast<int>(points.size()));
    shape->setPrecisePoints(precisePoints.empty() ? NULL :
			    precisePoints.data(),
			    static_cast<int>(points.size()));
    shape->geometryKind = geometryKind && geometryKind[0] ?
			  geometryKind : "line";
    prepared.release()->unrefNoDelete();
    return shape;
}

static int
publish_primitive_realization_line_set(
    SoBRLDatabaseSource *source,
    struct rt_primitive_lod_realization *realization,
    const char *sourceType,
    BObolSourceChildEffects *effects)
{
    if (!source || !realization || !realization->has_line_set)
	return 0;
    if (!realization->line_count) {
	primitive_realization_line_set_free(realization);
	return publish_external_primary_clear(*source, effects);
    }
    SoBRLVListShape *shape = vlist_from_primitive_realization_line_set(realization, "line");
    if (!shape)
	return -1;
    SbModernUtils::SoNodeRef prepared(shape);
    assign_external_primary_identity(shape, source,
	external_string_or_default(sourceType, "primitive-wireframe"), "line");
    const SoBRLVListShape *geometry = shape->getGeometrySource();
    return publish_prepared_external_shape(*source, *shape,
	external_bounds(geometry->point.getNum() ? geometry->point.getValues(0) : nullptr,
	    geometry->point.getNum()), effects);
}

static void
configure_auxiliary_line_set(SoBRLVListShape *shape,
    const SoBRLDatabaseSource *source, const char *name,
    const SbVec3f *points, const int32_t *commands, int count,
    const BObolAuxiliaryLineSetDisplayState *displayState)
{
    const char *sourcePath = source->path.getValue().getString();
    const uint32_t revision = source->sourceRevision.getValue();
    SbString identity = source_record_identity(source, sourcePath);
    if (identity.getLength() > 0)
	identity += "::";
    identity += name;

    shape->sourcePath = sourcePath ? sourcePath : "";
    shape->sourceName = name;
    shape->sourceType = "auxiliary-line-set";
    shape->sourceId = revision;
    shape->displayName = name;
    shape->geometryName = name;
    shape->sourceIdentity = identity;
    shape->cacheIdentity = record_identity_with_revision(identity.getString(),
			   revision);
    shape->databaseIntent = TRUE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = FALSE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = FALSE;
    shape->drawMode = source_record_draw_mode(source);
    shape->recordRole = "auxiliary";
    shape->geometryKind = "line";
    if (displayState && displayState->valid) {
	shape->drawMode = displayState->drawMode;
	shape->visible = displayState->visible;
	shape->highlighted = displayState->highlighted;
	shape->lineStyle = displayState->lineStyle;
	shape->lineWidth = displayState->lineWidth;
	shape->transparency = displayState->transparency;
	shape->materialColorValid = displayState->materialColorValid;
	shape->materialColor = displayState->materialColor;
	shape->materialRevision = displayState->materialRevision;
	shape->colorOverride = FALSE;
	shape->color = displayState->materialColor;
    } else {
	sync_shape_display_state(shape, source);
    }
    sync_shape_placement_state(shape, source);
    shape->setLineSet(points, commands, count);
    sync_shape_owner_state(shape, source);
}

struct primitive_submodel_publish_ctx {
    SoBRLDatabaseSource *source = nullptr;
    std::vector<SbModernUtils::SoNodeRef> children;
    std::unordered_map<std::string, size_t> occurrenceCounts;
    SbBox3f bounds;
    bool failed = false;
    std::exception_ptr failure;
};

static union tree *
primitive_submodel_wireframe_leaf(struct db_tree_state *tsp,
    const struct db_full_path *pathp, struct directory *dp, void *clientData)
{
    auto *ctx = static_cast<primitive_submodel_publish_ctx *>(clientData);
    if (!ctx || !ctx->source || ctx->failed || !tsp || !tsp->ts_dbip || !dp)
	return TREE_NULL;
    /* Keep exceptions inside the C walk so it can retire traversal resources. */
    try {
	std::unique_ptr<char, void (*)(char *)> pathName(
	    pathp && pathp->fp_len > 0 ? db_path_to_string(pathp) : nullptr,
	    [](char *name) { bu_free(name, "BObol submodel leaf path string"); });
	std::string name = pathName && pathName.get()[0] ? pathName.get() :
	    (dp->d_namep ? dp->d_namep : "submodel_leaf");
	/* A borrowed database may have combination instance IDs disabled. Use
	 * the same duplicate-ordinal suffix as compact occurrence publication. */
	const size_t ordinal = ctx->occurrenceCounts[name]++;
	if (ordinal) name += "@" + std::to_string(ordinal);
	struct rt_db_internal intern;
	RT_DB_INTERNAL_INIT(&intern);
	if (rt_db_get_internal(&intern, dp, tsp->ts_dbip, NULL) < 0) {
	    ctx->failed = true;
	    return TREE_NULL;
	}
	std::unique_ptr<rt_db_internal, decltype(&rt_db_free_internal)> imported(&intern, rt_db_free_internal);
	std::vector<SbVec3f> points;
	std::vector<int32_t> commands;
	if (!internal_payload_magic_valid(&intern) ||
	    plot_internal_to_vlist_geometry(points, commands, &intern, tsp->ts_ttol, tsp->ts_tol) < 0 ||
	    points.empty() || points.size() != commands.size() || points.size() > size_t(INT_MAX)) {
	    ctx->failed = true;
	    return TREE_NULL;
	}
	for (SbVec3f &point : points) {
	    point_t original, transformed;
	    VSET(original, point[X], point[Y], point[Z]);
	    MAT4X3PNT(transformed, tsp->ts_mat, original);
	    point.setValue(float(transformed[X]), float(transformed[Y]), float(transformed[Z]));
	    ctx->bounds.extendBy(point);
	}
	SoBRLVListShape *shape = new SoBRLVListShape;
	SbModernUtils::SoNodeRef prepared(shape);
	shape->setName(SbName(name.c_str()));
	configure_auxiliary_line_set(shape, ctx->source, name.c_str(), points.data(), commands.data(),
	    int(points.size()), nullptr);
	prepare_realized_shape_owner(*shape, *ctx->source);
	ctx->children.push_back(std::move(prepared));
	return make_nop_tree();
    } catch (...) {
	ctx->failed = true;
	ctx->failure = std::current_exception();
	return TREE_NULL;
    }
}

static int
publish_primitive_submodel_wireframe(SoBRLDatabaseSource *source,
    struct rt_db_internal *intern, const struct bg_tess_tol *ttol,
    const struct bn_tol *tol, BObolSourceChildEffects *effects)
{
    if (!source || !intern || intern->idb_type != ID_SUBMODEL || !intern->idb_ptr)
	return 0;
    auto *submodel = static_cast<struct rt_submodel_internal *>(intern->idb_ptr);
    RT_SUBMODEL_CK_MAGIC(submodel);
    struct db_i *dbip = DBI_NULL;
    std::unique_ptr<db_i, decltype(&db_close)> opened(nullptr, db_close);
    if (bu_vls_strlen(&submodel->file) != 0) {
	dbip = db_open(bu_vls_addr(&submodel->file), DB_OPEN_READONLY);
	if (dbip == DBI_NULL) return -1;
	opened.reset(dbip);
	if (!db_is_directory_non_empty(dbip) && db_dirbuild(dbip) < 0) return -1;
    } else {
	RT_CK_DBI(submodel->dbip);
	dbip = const_cast<struct db_i *>(submodel->dbip);
    }
    struct bn_tol localTol;
    if (!tol) {
	BN_TOL_INIT_SET_TOL(&localTol);
	tol = &localTol;
    }
    struct bg_tess_tol localTtol;
    if (!ttol) {
	BG_TESS_TOL_INIT_SET_TOL(&localTtol);
	ttol = &localTtol;
    }
    struct db_tree_state state;
    RT_DBTS_INIT(&state);
    state.ts_dbip = dbip;
    state.ts_ttol = ttol;
    state.ts_tol = tol;
    MAT_COPY(state.ts_mat, submodel->root2leaf);
    primitive_submodel_publish_ctx ctx;
    ctx.source = source;
    const char *argv[] = {bu_vls_addr(&submodel->treetop), nullptr};
    const int ret = db_walk_tree_leaf_instances(dbip, 1, argv, 1, &state, 0, nullptr,
	primitive_submodel_wireframe_leaf, &ctx);
    if (ctx.failure) std::rethrow_exception(ctx.failure);
    if (ret < 0 || ctx.failed) return -1;
    PreparedExternalPrimaryClear empty(*source);
    std::vector<SoNode *> children = empty.nodes;
    for (const auto &child : ctx.children) children.push_back(child.get());
    BObolPreparedSourcePublication publication(*source, children, ctx.bounds, false);
    empty.publish(publication, effects);
    return ctx.children.empty() ? 0 : 1;
}

int
SoBRLDatabaseSource::publishPrimitiveWireframe(
    struct rt_db_internal *intern,
    const struct bg_tess_tol *ttol,
    const struct bn_tol *tol)
{
    return this->publishPrimitiveWireframe(intern, ttol, tol, nullptr);
}

int
SoBRLDatabaseSource::publishPrimitiveWireframe(
    struct rt_db_internal *intern,
    const struct bg_tess_tol *ttol,
    const struct bn_tol *tol,
    BObolSourceChildEffects *effects)
{
    if (!intern)
	return 0;
    if (!internal_payload_magic_valid(intern))
	return -1;

    if (intern->idb_type == ID_SUBMODEL)
	return publish_primitive_submodel_wireframe(this, intern, ttol, tol, effects);

    if (intern->idb_meth && intern->idb_meth->ft_wireframe_line_set) {
	struct rt_primitive_lod_realization realization;
	memset(&realization, 0, sizeof(realization));
	std::unique_ptr<rt_primitive_lod_realization, decltype(&primitive_realization_line_set_free)> storage(
	    &realization, primitive_realization_line_set_free);
	struct bn_tol localTol;
	const struct bn_tol *useTol = tol;
	if (!useTol) {
	    BN_TOL_INIT_SET_TOL(&localTol);
	    useTol = &localTol;
	}
	struct bg_tess_tol localTtol;
	const struct bg_tess_tol *useTtol = ttol;
	if (!useTtol) {
	    BG_TESS_TOL_INIT_SET_TOL(&localTtol);
	    useTtol = &localTtol;
	}
	int ret = intern->idb_meth->ft_wireframe_line_set(&realization,
	    intern, useTtol, useTol);
	if (ret < 0 || !realization.has_line_set) {
	    return -1;
	}
	return publish_primitive_realization_line_set(this, storage.release(),
		"line-set", effects);
    }

    SoBRLVListShape *shape =
	vlist_from_plot_internal(intern, this, ttol, tol);
    if (!shape)
	return -1;

    SbModernUtils::SoNodeRef prepared(shape);
    const bool annotation = primitive_is_annotation(intern->idb_type,
			    primitive_type_label(intern));
    assign_external_primary_identity(shape, this,
				     annotation ? "annotation" : "line-set",
				     annotation ? "annotation" : "line");
    const SoBRLVListShape *geometry = shape->getGeometrySource();
    return publish_prepared_external_shape(*this, *shape,
	external_bounds(geometry->point.getNum() ? geometry->point.getValues(0) : nullptr,
	    geometry->point.getNum()), effects);
}

void
SoBRLDatabaseSource::syncCompactInstanceDisplayState(void)
{
    this->rebuildCompactInstanceDisplayState(TRUE);
}

Obol::InstanceStyle
compact_entry_style_from_source(const SoBRLDatabaseSource *source,
	const BObolCompactInstanceEntry &entry, SbBool selected,
	SbBool highlighted)
{
    Obol::InstanceStyle style;
    if (!source)
	return style;

    const BObolRealizedShapeSummary &summary = entry.shapeSummary;
    const bool inheritColor = source->materialColorValid.getValue() &&
	(source->materialPolicy.getValue() != SoBRLDatabaseSource::MATERIAL_DATABASE ||
	 !entry.semantic.materialColorValid);
    style.hasColorOverride = true;
    cad_shape_color(selected, source->selectedColor.getValue(), highlighted,
	source->highlightedColor.getValue(), summary.ghosted,
	source->ghostedColor.getValue(), source->colorOverride.getValue(),
	source->color.getValue(), inheritColor || entry.semantic.materialColorValid,
	inheritColor ? source->materialColor.getValue() : entry.semantic.materialColor, source->color.getValue(),
	source->transparency.getValue(), style.color);
    style.lineWidth = source->lineWidth.getValue() > 0 ?
	static_cast<float>(source->lineWidth.getValue()) : 1.0f;
    if (source->lineStyle.getValue() != 0)
	style.linePattern = 0xcf33u;
    return style;
}

bool
compact_set_material_styles(BObolCompactInstanceEntry &entry, const SoBRLDatabaseSource &source, bool colorsOnly)
{
    const auto nextStyle = [&](const Obol::InstanceStyle &current, SbBool selected, SbBool highlighted) {
	const auto appearance = compact_entry_style_from_source(&source, entry, selected, highlighted);
	if (!colorsOnly) return appearance;
	auto result = current;
	for (int axis = 0; axis < 3; ++axis) result.color[axis] = appearance.color[axis];
	return result;
    };
    const auto normal = nextStyle(entry.normalStyle, FALSE, FALSE);
    const auto selected = nextStyle(entry.selectedStyle, TRUE, FALSE);
    const auto highlighted = nextStyle(entry.highlightedStyle, FALSE, TRUE);
    const bool changed = !compact_style_equal(entry.normalStyle, normal) ||
	!compact_style_equal(entry.selectedStyle, selected) || !compact_style_equal(entry.highlightedStyle, highlighted);
    entry.normalStyle = normal;
    entry.selectedStyle = selected;
    entry.highlightedStyle = highlighted;
    return changed;
}

static bool
compact_sync_entry_from_source(BObolCompactInstanceEntry &entry,
	const SoBRLDatabaseSource *source)
{
    if (!source) return false;

    // Display policy must not overwrite the occurrence's database semantics.
    BObolRealizedShapeSummary &summary = entry.shapeSummary;
    if (summary.databaseIntent || !summary.nonDatabaseSource)
	summary.displayName = source->displayName.getValue().getLength() ? source->displayName.getValue() :
	    (summary.sourceName.getLength() ? summary.sourceName : summary.path);
    summary.drawMode = source_record_draw_mode(source);
    summary.hiddenLine = summary.drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    summary.materialRevision = source->materialRevision.getValue();
    summary.colorOverride = source->colorOverride.getValue();
    summary.color = source->color.getValue();
    return compact_set_material_styles(entry, *source);
}

static Obol::InstanceUpdate *
compact_instance_update(BObolCompactInstanceIndex &index, size_t ordinal)
{
    const auto &entry = index.entries[ordinal];
    if (ordinal < index.instances.size() && index.instances[ordinal].instance == entry.instance)
	return &index.instances[ordinal];
    auto found = std::find_if(index.instances.begin(), index.instances.end(),
	[&entry](const Obol::InstanceUpdate &update) { return update.instance == entry.instance; });
    return found != index.instances.end() ? &*found : nullptr;
}

void
compact_sync_instance_style(BObolCompactInstanceIndex &index, size_t ordinal)
{
    if (auto *update = compact_instance_update(index, ordinal))
	update->record.style = index.entries[ordinal].style;
}

void
SoBRLDatabaseSource::rebuildCompactInstanceDisplayState(
    SbBool syncSourceState)
{
    if (!this->d->compactIndex)
	return;

    this->d->compactIndex->hiddenInstances.clear();
    this->d->compactIndex->selectedInstances.clear();
    this->d->compactIndex->unpickableInstances.clear();

    for (size_t i = 0; i < this->d->compactIndex->entries.size(); i++) {
	BObolCompactInstanceEntry &entry = this->d->compactIndex->entries[i];
	const SbBool previousVisible = entry.visible;
	const SbBool previousSelectable = entry.selectable;
	const SbBool previousSelected = entry.selected;
	const SbBool previousHighlighted = entry.highlighted;
	const bool appearanceChanged = syncSourceState && compact_sync_entry_from_source(entry, this);
	const SbString &sourceInstanceKey = compact_instance_identity(entry);
	if (entry.semantic.sourceInstanceKey != sourceInstanceKey) {
	    entry.semantic.sourceInstanceKey = sourceInstanceKey;
	    compact_note_semantic_change(entry);
	}

	if (entry.visible != previousVisible ||
	    entry.selectable != previousSelectable)
	    entry.visibilityRevision = compact_next_revision(
		entry.visibilityRevision);
	if (entry.selected != previousSelected ||
	    entry.highlighted != previousHighlighted)
	    entry.selectionRevision = compact_next_revision(
		entry.selectionRevision);
	if (appearanceChanged)
	    entry.appearanceRevision = compact_next_revision(
		entry.appearanceRevision);
	entry.style = compact_effective_style(entry);
	compact_sync_shape_summary_state(entry);
	if (!entry.visible)
	    this->d->compactIndex->hiddenInstances.push_back(entry.instance);
	if (entry.selected)
	    this->d->compactIndex->selectedInstances.push_back(entry.instance);
	if (!entry.selectable)
	    this->d->compactIndex->unpickableInstances.push_back(entry.instance);

	compact_sync_instance_style(*this->d->compactIndex, i);
    }
}

void
SoBRLDatabaseSource::syncCompactInstancePlacementState(void)
{
    if (!this->d->compactIndex)
	return;

    auto &index = *this->d->compactIndex;
    for (size_t ordinal = 0; ordinal < index.entries.size(); ++ordinal) {
	auto &entry = index.entries[ordinal];
	SbMatrix localTransform = entry.geometryTransform;
	localTransform.multRight(entry.placementTransform);
	entry.localTransform = localTransform;
	const SbMatrix nextMatrix = cad_instance_matrix(this,
	    entry.localTransform);
	if (!entry.localToSource.equals(nextMatrix, 0.000001f)) {
	    entry.localToSource = nextMatrix;
	    entry.placementRevision = compact_next_revision(
		entry.placementRevision);
	}
	if (auto *update = compact_instance_update(index, ordinal))
	    update->record.localToRoot = entry.localToSource;
    }
}

void
SoBRLDatabaseSource::syncRealizedShapeOwnerState(void)
{
    this->syncCompactInstanceDisplayState();
    this->markCompiledAssemblyDirty();
    sync_realized_shape_owner_state_in_node(this, this);
}

SoBRLVListShape *
SoBRLDatabaseSource::getRealizedShape(void) const
{
    return this->getRealizedShape(0);
}

static SoBRLVListShape *
find_shape_in_node(SoNode *node, int &index)
{
    if (!node)
	return NULL;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	if (index == 0)
	    return static_cast<SoBRLVListShape *>(node);
	index--;
	return NULL;
    }

    if (node->isOfType(SoGroup::getClassTypeId()) &&
	!node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    SoBRLVListShape *shape = find_shape_in_node(group->getChild(i), index);
	    if (shape)
		return shape;
	}
    }

    return NULL;
}

SoBRLVListShape *
SoBRLDatabaseSource::getRealizedShape(int index) const
{
    if (index < 0)
	return NULL;
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoBRLVListShape *shape = find_shape_in_node(this->getChild(i), index);
	if (shape)
	    return shape;
    }
    return NULL;
}

static SbBool
vlist_shape_is_auxiliary(const SoBRLVListShape *shape)
{
    return shape &&
	   bu_strcmp(shape->recordRole.getValue().getString(), "auxiliary") == 0 ?
	   TRUE : FALSE;
}


SoBRLVListShape *
SoBRLDatabaseSource::findAuxiliaryVListShape(const char *name) const
{
    if (!name || !name[0])
	return NULL;

    for (int i = 0; i < this->getNumChildren(); i++) {
	SoNode *node = this->getChild(i);
	if (!node || !node->isOfType(SoBRLVListShape::getClassTypeId()))
	    continue;

	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	if (vlist_shape_is_auxiliary(shape) &&
	    bu_strcmp(shape->geometryName.getValue().getString(), name) == 0)
	    return shape;
    }

    return NULL;
}

SoBRLDatabaseSource *
SoBRLDatabaseSource::findAuxiliarySource(const char *sourcePath) const
{
    if (!sourcePath || !sourcePath[0])
	return NULL;

    for (int i = 0; i < this->getNumChildren(); i++) {
	SoNode *node = this->getChild(i);
	if (!node || !node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	    continue;

	SoBRLDatabaseSource *source = static_cast<SoBRLDatabaseSource *>(node);
	if (!source->auxiliarySource.getValue())
	    continue;

	const char *candidate = source->path.getValue().getString();
	if (bu_strcmp(candidate, sourcePath) == 0 ||
	    bu_strcmp(database_source_skip_leading_slash(candidate),
		   database_source_skip_leading_slash(sourcePath)) == 0)
	    return source;
    }

    return NULL;
}


static int
publish_auxiliary_removal(SoBRLDatabaseSource &source, std::vector<int> indices,
    BObolSourceChildEffects *effects)
{
    if (indices.empty()) return 0;
    const int removed = static_cast<int>(indices.size());
    auto removal = source.getChildren()->prepareRemoval(indices);
    std::vector<SoNode *> nextChildren;
    nextChildren.reserve(size_t(source.getNumChildren()) - indices.size());
    for (int i = 0; i < source.getNumChildren(); ++i)
	if (std::find(indices.begin(), indices.end(), i) == indices.end())
	    nextChildren.push_back(source.getChild(i));
    if (effects) {
	effects->stageChildOrder(source, nextChildren);
	effects->prepare();
    }
    removal->commit();
    if (effects) effects->commit();
    removal->notify();
    return removed;
}

static int
remove_auxiliary_child(SoBRLDatabaseSource &source, SoNode *child,
    BObolSourceChildEffects *effects)
{
    if (!child) return 0;
    const int index = source.findChild(child);
    if (index < 0) return 0;
    return publish_auxiliary_removal(source, {index}, effects);
}

static void
copy_auxiliary_shape_state(SoBRLVListShape &next, const SoBRLVListShape &previous)
{
    /* Copy semantic/display overrides, not geometry that this operation replaces.
     * In particular, primitive selection and edit intent survive named updates. */
    next.setName(previous.getName());
    const auto geometryFields = external_geometry_fields(next);
    const SoFieldData *data = static_cast<const SoFieldContainer &>(next).getFieldData();
    for (int i = 0; i < data->getNumFields(); ++i) {
	SoField *field = data->getField(&next, i);
	const SoField *old = previous.getField(data->getFieldName(i));
	if (field != &next.sharedGeometry &&
	    std::find(geometryFields.begin(), geometryFields.end(), field) == geometryFields.end())
	    field->copyFrom(*old);
	field->enableNotify(old->isNotifyEnabled());
    }
    next.enableNotify(previous.isNotifyEnabled());
}

class PreparedAuxiliaryLine {
public:
    PreparedAuxiliaryLine(SoBRLDatabaseSource &source, const SoBRLDatabaseSource &configuration,
	const char *name, const SbVec3f *points, const int32_t *commands, int count,
	const BObolAuxiliaryLineSetDisplayState *displayState) :
	previous(source.findAuxiliaryVListShape(name)), owner(new SoBRLVListShape)
    {
	auto &next = *static_cast<SoBRLVListShape *>(this->owner.get());
	next.setName(SbName(name));
	if (this->previous) copy_auxiliary_shape_state(next, *this->previous);
	std::vector<int32_t> fallbackCommands;
	if (!commands) {
	    fallbackCommands.assign(size_t(count), SoBRLVListShape::DRAW);
	    fallbackCommands.front() = SoBRLVListShape::MOVE;
	    commands = fallbackCommands.data();
	}
	configure_auxiliary_line_set(&next, &configuration, name, points, commands, count, displayState);
	std::vector<SoNode *> nextChildren;
	nextChildren.reserve(size_t(source.getNumChildren()) + (this->previous ? 0 : 1));
	for (int i = 0; i < source.getNumChildren(); ++i) {
	    SoNode *child = source.getChild(i);
	    nextChildren.push_back(child == this->previous ? &next : child);
	}
	if (!this->previous) nextChildren.push_back(&next);
	this->publication = std::make_unique<PreparedSourceChildren>(source, configuration, std::move(nextChildren));
    }
    SoBRLVListShape *replacedShape() const { return this->previous; }
    const std::vector<SoNode *> &order() const { return this->publication->order(); }
    void commit() { this->publication->commit(); }
    void restore() { this->publication->restore(); }
    void notify(std::exception_ptr &failure) { this->publication->notify(failure); }
    void notify() { this->publication->notify(); }
private:
    SoBRLVListShape *previous;
    SbModernUtils::SoNodeRef owner;
    std::unique_ptr<PreparedSourceChildren> publication;
};

int
SoBRLDatabaseSource::setAuxiliaryLineSet(const char *name,
    const SbVec3f *points, const int32_t *commands, int count,
    const BObolAuxiliaryLineSetDisplayState *displayState)
{
    return this->publishAuxiliaryLineSet(name, points, commands, count, displayState, nullptr);
}

int
SoBRLDatabaseSource::publishAuxiliaryLineSet(const char *name,
    const SbVec3f *points, const int32_t *commands, int count,
    const BObolAuxiliaryLineSetDisplayState *displayState,
    BObolSourceChildEffects *effects)
{
    if (!name || !name[0] || count < 0 || (count > 0 && !points))
	return 0;
    if (!count)
	return remove_auxiliary_child(*this, this->findAuxiliaryVListShape(name), effects);
    PreparedAuxiliaryLine publication(*this, *this, name, points, commands, count, displayState);
    if (effects) {
	effects->stageChildOrder(*this, publication.order());
	effects->prepare();
    }
    publication.commit();
    if (effects) effects->commit();
    publication.notify();
    return 1;
}


int
SoBRLDatabaseSource::setAuxiliarySourceLineSet(const char *sourcePath,
    const char *auxDisplayName, const SbVec3f *points, const int32_t *commands,
    int count, const BObolAuxiliaryLineSetDisplayState *displayState)
{
    return this->publishAuxiliarySourceLineSet(sourcePath, auxDisplayName, points,
	commands, count, displayState, nullptr);
}

int
SoBRLDatabaseSource::publishAuxiliarySourceLineSet(const char *sourcePath,
    const char *auxDisplayName, const SbVec3f *points, const int32_t *commands,
    int count, const BObolAuxiliaryLineSetDisplayState *displayState,
    BObolSourceChildEffects *effects)
{
    if (!sourcePath || !sourcePath[0] || count < 0 || (count > 0 && !points))
	return 0;
    SoBRLDatabaseSource *source = this->findAuxiliarySource(sourcePath);
    if (!count)
	return remove_auxiliary_child(*this, source, effects);

    const bool inserted = source == nullptr;
    SbModernUtils::SoNodeRef sourceOwner(source ? source : new SoBRLDatabaseSource);
    source = static_cast<SoBRLDatabaseSource *>(sourceOwner.get());
    const std::string sourceNameStorage = stable_name_from_path(sourcePath, 1);
    const char *sourceName = sourceNameStorage.empty() ? sourcePath : sourceNameStorage.c_str();
    if (inserted) source->setName(SbName(sourceName));

    auto *candidate = source->createDetachedRealizationTemplate();
    SbModernUtils::SoNodeRef candidateOwner(candidate);
    candidate->unref();
    auto &next = *candidate;
    /* Only configuration is copied. The named source retains its routing
     * identity and field sensors; no live field changes during preparation. */
    next.instanceKey = sourcePath;
    next.path = sourcePath;
    next.representationKey = sourcePath;
    next.representationMode = REPRESENTATION_DEFAULT;
    next.drawMode = this->drawMode.getValue() == SHADED ? SHADED : WIREFRAME;
    next.sourceRevision = this->sourceRevision.getValue();
    (void)next.d->setDatabaseBinding(this->d->dbip);
    next.auxiliarySource = TRUE;
    next.displayName = (auxDisplayName && auxDisplayName[0]) ?
			  auxDisplayName : sourceName;
    next.visible = this->visible.getValue();
    next.selected = this->selected.getValue();
    next.highlighted = this->highlighted.getValue();
    next.lineStyle = this->lineStyle.getValue();
    next.lineWidth = this->lineWidth.getValue();
    next.transparency = this->transparency.getValue();
    next.materialColorValid = this->materialColorValid.getValue();
    next.materialColor = this->materialColor.getValue();
    next.materialRevision = this->materialRevision.getValue();
    next.databaseMetadataValid = this->databaseMetadataValid.getValue();
    next.databaseRegionId = this->databaseRegionId.getValue();
    next.databaseAirCode = this->databaseAirCode.getValue();
    next.databaseMaterialId = this->databaseMaterialId.getValue();
    next.databaseLos = this->databaseLos.getValue();
    next.databaseMaterialColorValid =
	this->databaseMaterialColorValid.getValue();
    next.databaseMaterialColor = this->databaseMaterialColor.getValue();
    next.databaseMaterialShader =
	this->databaseMaterialShader.getValue();
    next.colorOverride = this->colorOverride.getValue();
    next.color = this->color.getValue();
    next.selectedColor = this->selectedColor.getValue();
    next.highlightedColor = this->highlightedColor.getValue();
    next.ghostedColor = this->ghostedColor.getValue();
    next.drawMatrixValid = this->drawMatrixValid.getValue();
    next.drawMatrix = this->drawMatrix.getValue();
    next.drawCenterValid = this->drawCenterValid.getValue();
    next.drawCenter = this->drawCenter.getValue();
    next.drawSizeValid = this->drawSizeValid.getValue();
    next.drawSize = this->drawSize.getValue();

    next.realizedRevision = next.sourceRevision.getValue();
    next.realizedSourceRevision = next.sourceRevision.getValue();
    next.realizedInputsRevision = next.inputsRevision.getValue();
    next.realizedViewRevision = next.viewRevision.getValue();
    next.realizationStatus = REALIZED;
    next.realizationDiagnostic = "";
    next.realizationIdentity = source_realization_identity(&next);
    next.stale = FALSE;
    next.staleReason = STALE_NONE;
    next.clearSourceBounds();

    const char *shapeName = (auxDisplayName && auxDisplayName[0]) ? auxDisplayName : sourceName;
    PreparedAuxiliaryLine line(*source, next, shapeName, points, commands, count, displayState);
    std::vector<std::unique_ptr<PreparedOwnedMetadata>> shapes;
    std::unordered_set<SoNode *> visited;
    for (int i = 0; i < source->getNumChildren(); ++i) {
	SoNode *child = source->getChild(i);
	if (child != line.replacedShape())
	    prepare_owned_shape_edit(child, next, next.sourceRevision.getValue(), shapes, visited);
    }
    PreparedScalarFields fields(*source, next);
    std::vector<SoNode *> nextChildren;
    nextChildren.reserve(size_t(this->getNumChildren()) + (inserted ? 1 : 0));
    for (int i = 0; i < this->getNumChildren(); ++i) nextChildren.push_back(this->getChild(i));
    if (inserted) nextChildren.push_back(source);
    PreparedSourceChildren parent(*this, std::move(nextChildren));
    const auto retiredStream = source->d->compactStagedSourceStream;

    const auto commitSourceState = [&] {
	(void)source->d->setDatabaseBinding(this->d->dbip);
	source->d->compactSourceProfile = BObolCompactSourceProfile();
	source->d->displayMeshLodContractRevisionValid = FALSE;
	source->d->compactStagedSourceStream.reset();
	bobol_identity_advance(source->d->compactPopulationEpoch);
	fields.commit();
    };
    /* The new source is still private, so expose its prepared identity to the
     * scene walk without exposing a partial publication to observers. */
    const bool preparedInsertedState = inserted && effects;
    if (preparedInsertedState) commitSourceState();

    if (effects) {
	effects->stageChildOrder(*source, line.order());
	if (inserted) effects->stageChildOrder(*this, parent.order());
	for (const auto &shape : shapes)
	    if (shape->changed()) effects->stageFrameEffect(*shape->target());
	effects->prepare();
    }

    if (!preparedInsertedState) commitSourceState();
    for (auto &shape : shapes) shape->commit();
    line.commit();
    source->markCompiledAssemblyDirty();
    source->markCadBatchDirty();
    parent.commit();
    if (effects) effects->commit();
    fields.restore();
    for (auto &shape : shapes) shape->restore();
    line.restore();
    parent.restore();
    std::exception_ptr failure;
    fields.notify(failure);
    for (auto &shape : shapes) shape->notify(failure);
    line.notify(failure);
    parent.notify(failure);
    if (failure) std::rethrow_exception(failure);
    return 1;
}


int
SoBRLDatabaseSource::clearAuxiliaryShapes(void)
{
    return this->removeAuxiliaryShapes(nullptr);
}

int
SoBRLDatabaseSource::removeAuxiliaryShapes(BObolSourceChildEffects *effects)
{
    std::vector<int> indices;
    for (int i = 0; i < this->getNumChildren(); ++i) {
	SoNode *child = this->getChild(i);
	if (node_is_auxiliary_source(child) || node_is_auxiliary_vlist(child))
	    indices.push_back(i);
    }
    return publish_auxiliary_removal(*this, std::move(indices), effects);
}

SoBRLMeshShape *
SoBRLDatabaseSource::getRealizedMesh(void) const
{
    return this->getRealizedMesh(0);
}

static SoBRLMeshShape *
find_mesh_in_node(SoNode *node, int &index)
{
    if (!node)
	return NULL;

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	if (index == 0)
	    return static_cast<SoBRLMeshShape *>(node);
	index--;
	return NULL;
    }

    if (node->isOfType(SoGroup::getClassTypeId()) &&
	!node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    SoBRLMeshShape *shape = find_mesh_in_node(group->getChild(i), index);
	    if (shape)
		return shape;
	}
    }

    return NULL;
}

SoBRLMeshShape *
SoBRLDatabaseSource::getRealizedMesh(int index) const
{
    if (index < 0)
	return NULL;
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoBRLMeshShape *shape = find_mesh_in_node(this->getChild(i), index);
	if (shape)
	    return shape;
    }
    return NULL;
}

SoBRLMaterialObject *
SoBRLDatabaseSource::getRealizedMaterialObject(void) const
{
    return this->getRealizedMaterialObject(0);
}

static SoBRLMaterialObject *
find_material_object_in_node(SoNode *node, int &index)
{
    if (!node)
	return NULL;

    if (node->isOfType(SoBRLMaterialObject::getClassTypeId())) {
	if (index == 0)
	    return static_cast<SoBRLMaterialObject *>(node);
	index--;
	return NULL;
    }

    if (node->isOfType(SoGroup::getClassTypeId()) &&
	!node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    SoBRLMaterialObject *object =
		find_material_object_in_node(group->getChild(i), index);
	    if (object)
		return object;
	}
    }

    return NULL;
}

SoBRLMaterialObject *
SoBRLDatabaseSource::getRealizedMaterialObject(int index) const
{
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoBRLMaterialObject *object =
	    find_material_object_in_node(this->getChild(i), index);
	if (object)
	    return object;
    }
    return NULL;
}

static int
count_shapes_in_node(SoNode *node)
{
    if (!node)
	return 0;

    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return 1;

    int ret = 0;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    ret += count_shapes_in_node(group->getChild(i));
    }

    return ret;
}

int
SoBRLDatabaseSource::getRealizedShapeCount(void) const
{
    int ret = 0;
    for (int i = 0; i < this->getNumChildren(); i++)
	ret += count_shapes_in_node(this->getChild(i));
    return ret;
}

SbBool
SoBRLDatabaseSource::hasRealizedWireGeometry(void) const
{
    /*
     * Compact indexes maintain these channel counts transactionally on every
     * append and in-place geometry replacement.  Scanning all occurrences
     * here made each LoD submission timer O(scene size); at 50k leaves this
     * single predicate dominated the otherwise idle owner thread.
     */
    if (this->d->compactIndexActive && this->d->compactIndex)
	return this->d->compactIndex->wireCount > 0 ? TRUE : FALSE;
    return this->getRealizedShapeCount() > 0 ? TRUE : FALSE;
}

static int
count_meshes_in_node(SoNode *node)
{
    if (!node)
	return 0;

    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return 1;

    int ret = 0;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    ret += count_meshes_in_node(group->getChild(i));
    }

    return ret;
}

int
SoBRLDatabaseSource::getRealizedMeshCount(void) const
{
    int ret = 0;
    for (int i = 0; i < this->getNumChildren(); i++)
	ret += count_meshes_in_node(this->getChild(i));
    return ret;
}

SbBool
SoBRLDatabaseSource::hasRealizedMeshGeometry(void) const
{
    if (this->d->compactIndexActive && this->d->compactIndex)
	return this->d->compactIndex->shadedCount > 0 ? TRUE : FALSE;
    return this->getRealizedMeshCount() > 0 ? TRUE : FALSE;
}

SbBool
SoBRLDatabaseSource::hasDisplayMeshLodRequests(void) const
{
    return (this->d->compactIndexActive && this->d->compactIndex &&
	    this->d->compactIndex->sourceMeshRequestCount > 0 &&
	    this->d->displayMeshLodContractRevisionValid &&
	    this->d->displayMeshLodContractSourceRevision ==
		this->sourceRevision.getValue() &&
	    this->d->displayMeshLodContractInputsRevision ==
	this->inputsRevision.getValue()) ? TRUE : FALSE;
}

SbBool
SoBRLDatabaseSource::hasDisplayLodTargets(void) const
{
    return this->getDisplayLodTargetCount() > 0 ? TRUE : FALSE;
}

size_t
SoBRLDatabaseSource::getDisplayLodTargetCount(void) const
{
    if (!this->d->compactIndexActive || !this->d->compactIndex)
	return 0;
    if (this->hasDisplayMeshLodRequests())
	return this->d->compactIndex->displayLodTargetCount;
    return this->d->compactIndex->residentProgressiveGeometryCount;
}

SbBool
SoBRLDatabaseSource::hasDisplayResidentProgressiveGeometry(void) const
{
    return this->getDisplayResidentProgressiveGeometryCount() > 0 ?
	TRUE : FALSE;
}

size_t
SoBRLDatabaseSource::getDisplayResidentProgressiveGeometryCount(void) const
{
    return this->d->compactIndexActive && this->d->compactIndex ?
	this->d->compactIndex->residentProgressiveGeometryCount : 0;
}

size_t
SoBRLDatabaseSource::getDisplayMeshLodRequestCount(void) const
{
    return this->hasDisplayMeshLodRequests() ?
	this->d->compactIndex->sourceMeshRequestCount : 0;
}

uint64_t
SoBRLDatabaseSource::getDisplayMeshLodRevision(void) const
{
    return this->d->displayMeshLodRevision;
}

template <typename Delta>
static SbBool
database_source_display_mesh_lod_changed_entries(
    uint64_t revision, uint64_t currentRevision, uint64_t floorRevision,
    const std::deque<Delta> &deltas, std::vector<size_t> &entryIndices,
    SbBool *coverageInvalidated)
{
    entryIndices.clear();
    if (coverageInvalidated)
	*coverageInvalidated = FALSE;
    if (revision == currentRevision)
	return TRUE;
    if (!revision || revision > currentRevision || revision < floorRevision)
	return FALSE;

    for (const Delta &delta : deltas) {
	if (delta.revision <= revision)
	    continue;
	entryIndices.insert(entryIndices.end(), delta.entryIndices.begin(),
	    delta.entryIndices.end());
	if (coverageInvalidated && delta.coverageInvalidated)
	    *coverageInvalidated = TRUE;
    }
    if (entryIndices.empty())
	return FALSE;
    std::sort(entryIndices.begin(), entryIndices.end());
    entryIndices.erase(std::unique(entryIndices.begin(),
	entryIndices.end()), entryIndices.end());
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getDisplayMeshLodChangedEntries(
    uint64_t revision, std::vector<size_t> &entryIndices,
    SbBool *coverageInvalidated) const
{
    return database_source_display_mesh_lod_changed_entries(
	revision, this->d->displayMeshLodRevision,
	this->d->displayMeshLodDeltaFloorRevision,
	this->d->displayMeshLodDeltas, entryIndices, coverageInvalidated);
}

uint64_t
SoBRLDatabaseSource::getDisplayMeshLodVisibilityRevision(void) const
{
    return this->d->displayMeshLodVisibilityRevision;
}

SbBool
SoBRLDatabaseSource::getDisplayMeshLodVisibilityChangedEntries(
    uint64_t revision, std::vector<size_t> &entryIndices) const
{
    return database_source_display_mesh_lod_changed_entries(
	revision, this->d->displayMeshLodVisibilityRevision,
	this->d->displayMeshLodVisibilityDeltaFloorRevision,
	this->d->displayMeshLodVisibilityDeltas, entryIndices, NULL);
}


int
SoBRLDatabaseSource::refreshCompactObjectGeometry(
    const char *objectPath, uint32_t nextSourceRevision)
{
    if (!this->isCompactOccurrenceRegistry() || !this->d->compactIndex ||
	!this->d->dbip || !objectPath || !objectPath[0])
	return 0;

    const std::string databasePath =
	database_source_db_path_without_instance_suffixes(objectPath);
    const char *objectName = strrchr(databasePath.c_str(), '/');
    objectName = objectName && objectName[1] ? objectName + 1 : databasePath.c_str();
    if (!objectName[0])
	return 0;

    struct directory *dp = db_lookup(this->d->dbip, objectName, LOOKUP_QUIET);
    if (!dp)
	return -1;
    if (dp->d_flags & RT_DIR_COMB) {
	const uint32_t revision = nextSourceRevision ? nextSourceRevision :
	    bobol_identity_successor_or_terminate(this->sourceRevision.getValue());
	return BObolPreparedSourcePublication::refreshCombination(*this, revision);
    }

    std::vector<size_t> matching;
    const std::unordered_map<std::string, std::vector<size_t>>::const_iterator
	bySourceName = this->d->compactIndex->entryIndicesBySourceName.find(objectName);
    if (bySourceName != this->d->compactIndex->entryIndicesBySourceName.end())
	matching.insert(matching.end(), bySourceName->second.begin(),
	    bySourceName->second.end());
    const std::unordered_map<std::string, std::vector<size_t>>::const_iterator
	byLeaf = this->d->compactIndex->entryIndicesByLeaf.find(objectName);
    if (byLeaf != this->d->compactIndex->entryIndicesByLeaf.end())
	matching.insert(matching.end(), byLeaf->second.begin(), byLeaf->second.end());
    std::sort(matching.begin(), matching.end());
    matching.erase(std::unique(matching.begin(), matching.end()),
	matching.end());
    if (matching.empty())
	return 0;

    owned_leaf_internal validInternal;
    if (rt_db_get_internal(&validInternal.local, dp, this->d->dbip, NULL) < 0 ||
	!internal_payload_magic_valid(&validInternal.local)) {
	if (validInternal.local.idb_ptr)
	    rt_db_free_internal(&validInternal.local);
	return -1;
    }
    validInternal.ownsLocal = true;

    const char *typeLabel = primitive_type_label(&validInternal.local);
    const int sourceDrawMode = source_record_draw_mode(this);
    const bool wantWire = sourceDrawMode == BOBOL_LOD_DRAW_WIRE;
    SoBRLVListShape *sharedWire = NULL;
    SoBRLMeshShape *sharedMesh = NULL;
    Obol::PartGeometryBuilder generated;
    bool geometryValid = false;
    bool directWire = false;
    bool directMesh = false;
    bool replacementViewDependentCsgGeometry = false;
    bool replacementLodBacked = false;
    bool replacementSourceMeshRequestValid = false;
    BObolSourceMeshRequest replacementSourceMeshRequest;
    SbBox3f localBounds;
    (void)local_bounds_from_internal(&validInternal.local, localBounds);
    if (wantWire) {
	if (validInternal.local.idb_type == ID_BOT)
	    geometryValid = cad_wire_part_geometry_from_bot(
		static_cast<const struct rt_bot_internal *>(
		    validInternal.local.idb_ptr), generated) != 0;
	else
	    geometryValid = cad_wire_part_geometry_from_lod_realization_internal(
		&validInternal.local, this, localBounds, generated,
		&replacementViewDependentCsgGeometry) != 0;
	if (!geometryValid)
	    geometryValid = cad_wire_part_geometry_from_plot_internal(
		&validInternal.local, this, generated) != 0;
	directWire = geometryValid;
	if (!geometryValid && validInternal.local.idb_type == ID_BOT)
	    sharedWire = vlist_from_bot_wireframe(
		static_cast<const struct rt_bot_internal *>(
		    validInternal.local.idb_ptr));
	else if (!geometryValid)
	    sharedWire = vlist_from_lod_realization_internal(
		&validInternal.local, this, localBounds);
	if (!geometryValid && !sharedWire)
	    sharedWire = vlist_from_plot_internal(&validInternal.local, this);
    } else {
	const struct rt_bot_internal *bot =
	    validInternal.local.idb_type == ID_BOT ?
	    static_cast<const struct rt_bot_internal *>(
		validInternal.local.idb_ptr) : NULL;
	if (bot && this->lodBotThreshold.getValue() > 0 &&
	    bot->num_faces >= this->lodBotThreshold.getValue() &&
	    cad_source_mesh_request_from_bot(replacementSourceMeshRequest,
		bot) &&
	    cad_wire_part_geometry_from_aabb(
		replacementSourceMeshRequest.bounds, generated)) {
	    geometryValid = true;
	    replacementLodBacked = true;
	    replacementSourceMeshRequestValid = true;
	} else {
	    geometryValid = cad_mesh_part_geometry_from_internal(
		&validInternal.local, this, generated) != 0;
	}
	if (geometryValid &&
	    sourceDrawMode == BOBOL_LOD_DRAW_HIDDEN_LINE)
	    (void)cad_mesh_append_hidden_line_edges(generated);
	directMesh = geometryValid;
	if (!geometryValid)
	    sharedMesh = mesh_from_internal(&validInternal.local, this);
	if (!geometryValid && !sharedMesh) {
	    sharedWire = vlist_from_lod_realization_internal(
		&validInternal.local, this, localBounds);
	    if (!sharedWire)
		sharedWire = vlist_from_plot_internal(&validInternal.local, this);
	}
    }
    if (!geometryValid && !sharedWire && !sharedMesh)
	return -1;

    SbModernUtils::SoNodeRef wireOwner(sharedWire), meshOwner(sharedMesh);
    const uint32_t revision = nextSourceRevision ? nextSourceRevision :
	bobol_identity_successor_or_terminate(
	    this->sourceRevision.getValue());
    if (sharedWire)
	assign_shared_geometry_identity(sharedWire, objectName, typeLabel,
	    revision, "line");
    if (sharedMesh)
	assign_shared_geometry_identity(sharedMesh, objectName, typeLabel,
	    revision, "surface");

    for (size_t index : matching) {
	const BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[index];
	if (((sharedWire || directWire) && !entry.wireGeometry &&
	     !entry.pointGeometry) ||
	    ((sharedMesh || directMesh) && !entry.meshGeometry &&
	     !entry.sourceMeshRequestValid)) {
	    return -1;
	}
    }
    if (!geometryValid && sharedWire) {
	SbModernUtils::SoNodeRef shapeOwner(new SoBRLVListShape);
	auto *geometryShape = static_cast<SoBRLVListShape *>(shapeOwner.get());
	geometryShape->setSharedGeometry(sharedWire);
	geometryValid = cad_vlist_part_geometry(geometryShape, generated) != 0;
    } else if (!geometryValid && sharedMesh) {
	SbModernUtils::SoNodeRef shapeOwner(new SoBRLMeshShape);
	auto *geometryShape = static_cast<SoBRLMeshShape *>(shapeOwner.get());
	geometryShape->setSharedGeometry(sharedMesh);
	geometryShape->hiddenLine =
	    sourceDrawMode == BOBOL_LOD_DRAW_HIDDEN_LINE ? TRUE : FALSE;
	geometryShape->drawMode = sourceDrawMode;
	geometryValid = cad_mesh_part_geometry(geometryShape, generated) != 0;
    }
    if (!geometryValid) {
	return -1;
    }

    std::string partKey;
    const char *partKind = compact_geometry_kind(generated);
    if (!cad_part_key_for_geometry(partKind,
	generated, partKey)) {
	return -1;
    }
    BObolCompactInstanceIndex &registry = *this->d->compactIndex;
    std::shared_ptr<const Obol::PartGeometry> geometry;
    const auto existingPart = registry.partIdByKey.find(partKey);
    if (existingPart != registry.partIdByKey.end()) {
	for (const auto &part : registry.parts)
	    if (part.part == existingPart->second) {
		geometry = part.geometry;
		break;
	    }
    }
    if (!geometry)
	geometry = bobol_cad_build_geometry(std::move(generated), "compact occurrence geometry");
    if (!geometry)
	return -1;
    CompactPartPublication preparedPart(registry, geometry);
    if (!preparedPart.valid())
	return -1;
    const Obol::PartId partId = preparedPart.part();
    CompactPublicationSlot newReference(registry.partReferenceCounts, partId);

    std::vector<BObolCompactInstanceEntry> replacements;
    replacements.reserve(matching.size());
    std::unordered_map<Obol::PartId, size_t, std::hash<Obol::PartId>> replacedPartCounts;
    bool rebuildBounds = registry.sourceBoundsDirty;
    for (size_t index : matching) {
	if (index >= registry.entries.size() || index >= registry.instances.size())
	    throw std::logic_error("compact edit has no matching instance record");
	const BObolCompactInstanceEntry &previous = registry.entries[index];
	++replacedPartCounts[previous.part];
	rebuildBounds = rebuildBounds ||
	    compact_bounds_removal_requires_rebuild(registry.sourceBounds, previous.sourceBounds);
	BObolCompactInstanceEntry entry = previous;
	/* The replacement is native object-local geometry.  Preserve the
	 * explicitly retained tree placement while discarding the old proxy/PCA
	 * geometry transform.  Inferring this with an inverse is invalid for
	 * degenerate boxes and needlessly couples replacement to payload shape. */
	entry.geometryTransform = SbMatrix::identity();
	entry.localTransform = entry.placementTransform;
	entry.localToSource = cad_instance_matrix(this, entry.localTransform);
	entry.part = partId;
	entry.geometry = geometry;
	entry.wireGeometry = geometry->wire ? TRUE : FALSE;
	entry.pointGeometry = geometry->points ? TRUE : FALSE;
	entry.meshGeometry = geometry->shaded ? TRUE : FALSE;
	entry.viewDependentCsgGeometry =
	    replacementViewDependentCsgGeometry ? TRUE : FALSE;
	entry.lodBacked = replacementLodBacked ? TRUE : FALSE;
	entry.shapeSummary.geometryKind =
	    replacementLodBacked ? "aabb" :
	    ((sharedWire || directWire) ? "line" : "surface");
	entry.sourceMeshRequestValid = replacementSourceMeshRequestValid ?
	    TRUE : FALSE;
	if (entry.sourceMeshRequestValid)
	    entry.sourceMeshRequest = replacementSourceMeshRequest;
	else
	    entry.sourceMeshRequest.clear();
	entry.geometryRevision = compact_next_revision(entry.geometryRevision);
	entry.semantic.sourceId = revision;
	compact_note_semantic_change(entry);
	compact_sync_shape_summary(entry);
	entry.shapeSummary.ownerSourceRevision = revision;
	entry.shapeSummary.cacheIdentity = record_identity_with_revision(
	    entry.shapeSummary.sourceIdentity.getString(), revision);
	if (entry.sourceMeshRequestValid) {
	    compact_source_mesh_request_sync(entry.sourceMeshRequest,
		entry.shapeSummary);
	    compact_summary_lod_from_source_mesh_request(entry.shapeSummary,
		entry.sourceMeshRequest);
	}
	entry.sourceBounds = database_source_transform_bounds(
	    compact_part_geometry_bounds(entry.geometry), entry.localTransform);
	replacements.push_back(std::move(entry));
    }

    std::unordered_set<Obol::PartId, std::hash<Obol::PartId>> releasedParts;
    for (const auto &replaced : replacedPartCounts) {
	const auto count = registry.partReferenceCounts.find(replaced.first);
	if (count == registry.partReferenceCounts.end() || count->second < replaced.second)
	    throw std::logic_error("compact edit has inconsistent part references");
	if (count->second == replaced.second && replaced.first != partId)
	    releasedParts.insert(replaced.first);
    }
    SbBox3f bounds = rebuildBounds ? SbBox3f() : registry.sourceBounds;
    if (rebuildBounds) {
	size_t replacement = 0;
	for (size_t i = 0; i < registry.entries.size(); ++i) {
	    const auto &entry = replacement < matching.size() && matching[replacement] == i ?
		replacements[replacement++] : registry.entries[i];
	    if (!entry.sourceBounds.isEmpty()) bounds.extendBy(entry.sourceBounds);
	}
    } else {
	for (const auto &entry : replacements)
	    if (!entry.sourceBounds.isEmpty()) bounds.extendBy(entry.sourceBounds);
    }
    PreparedSourceRealization publication(*this, revision, bounds, TRUE, false,
	this->realizationRoleFlags.getValue());
    std::vector<std::unique_ptr<PreparedOwnedMetadata>> shapes;
    std::unordered_set<SoNode *> visited;
    for (int i = 0; i < this->getNumChildren(); ++i)
	prepare_owned_shape_edit(this->getChild(i), *this, revision, shapes, visited);

    /* Commit only after entry, metadata and table preparation is complete.
     * Compiled presentation remains retained for its existing incremental
     * refresh; its active flag and all source facts change before callbacks. */
    preparedPart.commit();
    for (size_t i = 0; i < matching.size(); ++i)
	compact_commit_entry_geometry(registry, matching[i], std::move(replacements[i]));
    newReference.commit();
    for (const auto &released : releasedParts)
	registry.partReferenceCounts.erase(released);
    registry.parts.erase(
	std::remove_if(registry.parts.begin(),
	    registry.parts.end(),
	    [&](const BObolCompactPartReference &part) {
		return releasedParts.find(part.part) != releasedParts.end();
	    }), registry.parts.end());
    for (auto it = registry.partIdByKey.begin();
	 it != registry.partIdByKey.end();) {
	if (releasedParts.find(it->second) != releasedParts.end())
	    it = registry.partIdByKey.erase(it);
	else
	    ++it;
    }
    for (auto it = registry.partIdByGeometry.begin();
	 it != registry.partIdByGeometry.end();) {
	if (releasedParts.find(it->second.part) != releasedParts.end())
	    it = registry.partIdByGeometry.erase(it);
	else
	    ++it;
    }

    /* Encoded size and asset reuse were certified before this database edit. */
    this->d->compactSourceProfile = BObolCompactSourceProfile();
    registry.sourceBounds = bounds;
    registry.sourceBoundsDirty = false;
    this->d->displayMeshLodContractRevisionValid = registry.sourceMeshRequestCount > 0;
    this->d->displayMeshLodContractSourceRevision = revision;
    this->d->displayMeshLodContractInputsRevision = this->inputsRevision.getValue();
    publication.commit();
    for (auto &shape : shapes) shape->commit();
    this->markCompiledAssemblyDirty();
    this->markCadBatchDirty();
    this->markDisplayMeshLodDirty();
    publication.restore();
    for (auto &shape : shapes) shape->restore();
    std::exception_ptr failure;
    publication.notify(failure);
    for (auto &shape : shapes) shape->notify(failure);
    if (failure) std::rethrow_exception(failure);
    return static_cast<int>(matching.size());
}
