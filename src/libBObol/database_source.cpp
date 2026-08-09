/*                D A T A B A S E _ S O U R C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

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
#include "BObol/BVListShape.h"
#include "cad_assembly_private.h"
#include "compact_occurrence_registry_private.h"
#include "database_source_private.h"
#include "database_source_realization.h"
#include "performance_private.h"

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
#include "bu/parallel.h"
#include "bu/str.h"
#include "bu/time.h"
#include "bu/vls.h"
#include "nmg.h"
#include "raytrace.h"
#include "rt/func.h"
#include "rt/global.h"
#include "rt/db4.h"
#include "rt/nongeom.h"
#include "rt/db_fullpath.h"
#include "rt/eval_wireframe.h"
#include "rt/primitives/annot.h"
#include "rt/tree.h"
#include "rt/vlist.h"
#include "rt/view.h"
#include "wdb.h"

#include <Inventor/SbName.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/nodes/SoGroup.h>

#include <Inventor/nodes/SoMatrixTransform.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/sensors/SoFieldSensor.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <condition_variable>
#include <deque>
#include <inttypes.h>
#include <limits.h>
#include <limits>
#include <map>
#include <math.h>
#include <memory>
#include <mutex>
#include <numeric>
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

static uint64_t database_source_handle_id(void);

struct SoBRLDatabaseSource::Impl :
    BObolCadSourceState,
    BObolCompactOccurrenceRegistryState,
    BObolCadPresentationBridgeState,
    BObolDatabaseSourceSensorState
{
    Impl(void) :
	BObolCadSourceState(
	    database_source_handle_id(), database_source_handle_id())
    {
    }
};

static void realized_vlist_shape_summary(const SoBRLVListShape *shape,
	BObolRealizedShapeSummary &summary);
static void realized_mesh_shape_summary(const SoBRLMeshShape *shape,
	BObolRealizedShapeSummary &summary);
static int cad_source_mesh_request_from_bot(
	BObolSourceMeshRequest &request, const struct rt_bot_internal *bot);

static const SbString &
compact_instance_identity(const BObolCompactInstanceEntry &entry)
{
    if (entry.instanceKey.getLength() == 0)
	entry.instanceKey.sprintf("obol-instance:%016" PRIx64 "%016" PRIx64,
	    entry.instance.w1, entry.instance.w0);
    return entry.instanceKey;
}

static SbBox3f
compact_part_geometry_bounds(
    const std::shared_ptr<const Obol::PartGeometry> &geometry)
{
    SbBox3f bounds;
    bounds.makeEmpty();
    if (!geometry)
	return bounds;
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

static uint64_t
database_source_handle_id(void)
{
    uint64_t id = database_source_next_handle_id.fetch_add(1);
    if (id == 0)
	id = database_source_next_handle_id.fetch_add(1);
    return id;
}

static int
database_source_float_different(float a, float b)
{
    return fabsf(a - b) > 1.0e-6f;
}

static int
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

static int
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

static SbBox3f
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

static std::string
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
bobol_database_source_fullpath_material_color(
    struct db_i *dbip,
    const struct db_full_path *pathp,
    SbColor &color)
{
    if (!dbip || !pathp || pathp->fp_len <= 0)
	return FALSE;

    struct bu_color dbColor = BU_COLOR_INIT_ZERO;
    db_full_path_color(&dbColor, const_cast<struct db_full_path *>(pathp),
		       dbip);
    unsigned char rgb[3] = {255, 255, 255};
    bu_color_to_rgb_chars(&dbColor, rgb);
    color = SbColor(static_cast<float>(rgb[0]) / 255.0f,
		    static_cast<float>(rgb[1]) / 255.0f,
		    static_cast<float>(rgb[2]) / 255.0f);
    return TRUE;
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

struct BObolMaterialCombState {
    BObolMaterialCombState(void) :
	isCombination(false),
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
		this->applyCombination(state, dp);
		this->pathStates.emplace(prefix, state);
	    }

	    if (!slash)
		break;
	    cp = slash + 1;
	}

	if (!state.explicitColorValid && state.regionId >= 0)
	    state.color = this->regionColor(state.regionId);
	resolved = state;
	return true;
    }

private:
    const BObolMaterialCombState &combinationState(struct directory *dp)
    {
	std::unordered_map<struct directory *,
	    BObolMaterialCombState>::iterator found =
	    this->combStates.find(dp);
	if (found != this->combStates.end())
	    return found->second;

	BObolMaterialCombState state;
	state.isCombination = (dp->d_flags & RT_DIR_COMB) != 0;
	if (state.isCombination) {
	    struct rt_db_internal intern;
	    RT_DB_INTERNAL_INIT(&intern);
	    if (rt_db_get_internal(&intern, dp, this->dbip, NULL) >= 0) {
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
		rt_db_free_internal(&intern);
	    }
	}
	return this->combStates.emplace(dp, state).first->second;
    }

    void applyCombination(BObolMaterialPathState &pathState,
	struct directory *dp)
    {
	const BObolMaterialCombState &comb = this->combinationState(dp);
	if (!comb.isCombination)
	    return;

	if (!pathState.inRegion && pathState.colorInherit == DB_INH_LOWER &&
	    comb.rgbValid) {
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
    std::unordered_map<struct directory *, BObolMaterialCombState>
	combStates;
    std::unordered_map<std::string, BObolMaterialPathState> pathStates;
    std::unordered_map<int, SbColor> regionColors;
};

} // namespace


int
bobol_database_sources_refresh_material_colors(
    SoBRLDatabaseSource *const *sources,
    size_t sourceCount,
    uint32_t materialRevision,
    struct db_i *dbip)
{
    if (!sourceCount)
	return 0;
    if (!sources || !dbip)
	return -1;

    BObolMaterialColorSweep sweep(dbip);
    int changed = 0;
    for (size_t i = 0; i < sourceCount; i++) {
	SoBRLDatabaseSource *source = sources[i];
	if (!source || source->materialRevision.getValue() == materialRevision)
	    continue;
	if (source->hasCompactInstanceIndex()) {
	    if (source->refreshMaterialColorFromDatabase(materialRevision,
		    dbip) > 0)
		changed = 1;
	    continue;
	}

	BObolMaterialPathState resolved;
	if (!sweep.resolve(source->path.getValue().getString(), resolved))
	    continue;
	if (resolved.inRegion) {
	    (void)source->setDatabaseMetadataState(
		TRUE, resolved.regionId, resolved.airCode, resolved.materialId,
		resolved.los, TRUE, resolved.color,
		SbString(resolved.shader.c_str()));
	}
	const int sourceChanged = source->setDisplayState(
	    FALSE,
	    source->sourceRevision.getValue(),
	    source->inputsRevision.getValue(),
	    source->visible.getValue(),
	    source->selected.getValue(),
	    source->highlighted.getValue(),
	    source->lineStyle.getValue(),
	    source->lineWidth.getValue(),
	    source->transparency.getValue(),
	    source->colorOverride.getValue(),
	    source->color.getValue(),
	    TRUE,
	    resolved.color,
	    materialRevision);
	if (sourceChanged > 0)
	    changed = 1;
    }
    return changed;
}

BObolDatabaseSourceSummary::BObolDatabaseSourceSummary(void) :
    valid(FALSE),
    path(""),
    instanceKey(""),
    parentInstanceKey(""),
    occurrenceIndex(0),
    booleanOperation(SoBRLDatabaseSource::BOOLEAN_UNION),
    displayName(""),
    hasParent(FALSE),
    drawTreeDepth(0),
    parentGroupPath(""),
    representationKey(""),
    representationMode(-1),
    drawMode(SoBRLDatabaseSource::WIREFRAME),
    sourceRevision(0),
    inputsRevision(0),
    viewRevision(0),
    realizedRevision(0),
    realizedSourceRevision(0),
    realizedInputsRevision(0),
    realizedViewRevision(0),
    realizationStatus(SoBRLDatabaseSource::UNREALIZED),
    realizationDiagnostic(""),
    realizationIdentity(""),
    realizationRoleFlags(SoBRLDatabaseSource::REALIZATION_ROLE_NONE),
    realizationViewDependent(FALSE),
    realizationCsgLodEnabled(FALSE),
    realizationMeshLodEnabled(FALSE),
    realizationViewScale(0.0f),
    realizationLodScale(1.0f),
    realizationViewWidth(0),
    realizationViewHeight(0),
    realizationBotThreshold(0),
    realizationCurveScale(0.0f),
    realizationPointScale(0.0f),
    visible(TRUE),
    selected(FALSE),
    highlighted(FALSE),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0f),
    materialColorValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialRevision(0),
    materialPolicy(SoBRLDatabaseSource::MATERIAL_INHERIT),
    databaseMetadataValid(FALSE),
    databaseRegionId(0),
    databaseAirCode(0),
    databaseMaterialId(0),
    databaseLos(0),
    databaseMaterialColorValid(FALSE),
    databaseMaterialColor(1.0f, 1.0f, 1.0f),
    databaseMaterialShader(""),
    colorOverride(FALSE),
    color(1.0f, 1.0f, 1.0f),
    selectedColor(1.0f, 1.0f, 1.0f),
    highlightedColor(1.0f, 1.0f, 0.0f),
    ghostedColor(0.55f, 0.55f, 0.55f),
    drawMatrixValid(FALSE),
    drawMatrix(SbMatrix::identity()),
    drawCenterValid(FALSE),
    drawCenter(0.0f, 0.0f, 0.0f),
    drawSizeValid(FALSE),
    drawSize(0.0f),
    sourceBoundsValid(FALSE),
    sourceBounds(),
    stale(TRUE),
    staleReason(SoBRLDatabaseSource::STALE_SOURCE),
    realizedShapeCount(0),
    realizedMeshCount(0),
    realizedMaterialObjectCount(0)
{
    this->sourceBounds.makeEmpty();
}

BObolAuxiliaryLineSetDisplayState::BObolAuxiliaryLineSetDisplayState(void) :
    valid(FALSE),
    drawMode(BOBOL_LOD_DRAW_WIRE),
    visible(TRUE),
    highlighted(FALSE),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0f),
    materialColorValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialRevision(0)
{
}

BObolCompactInstanceHandle::BObolCompactInstanceHandle(void) :
    sourceNodeId(0),
    instanceWord0(0),
    instanceWord1(0)
{
}

SbBool
BObolCompactInstanceHandle::isValid(void) const
{
    return sourceNodeId != 0 && (instanceWord0 != 0 || instanceWord1 != 0) ?
	TRUE : FALSE;
}

BObolCompactInstanceSummary::BObolCompactInstanceSummary(void) :
    valid(FALSE),
    sourceContentHash(0),
    sourceFaceCount(0),
    sourcePointCount(0),
    localToSource(SbMatrix::identity()),
    geometryIdentity(0),
    geometryRevision(0),
    appearanceRevision(0),
    placementRevision(0),
    visibilityRevision(0),
    selectionRevision(0),
    occurrenceIndex(0),
    booleanOperation(SoBRLDatabaseSource::BOOLEAN_UNION),
    regionId(0),
    airCode(0),
    materialId(0),
    los(0),
    materialColorValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialShader(""),
    appearanceColorValid(FALSE),
    appearanceColor(1.0f, 1.0f, 1.0f),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0f),
    wireGeometry(FALSE),
    pointGeometry(FALSE),
    meshGeometry(FALSE),
    lodBacked(FALSE),
    sourceMeshRequestValid(FALSE),
    visible(TRUE),
    selectable(TRUE),
    selected(FALSE),
    highlighted(FALSE)
{
    meshAssetBounds.makeEmpty();
}

BObolCompactLodInstanceSummary::BObolCompactLodInstanceSummary(void) :
    valid(FALSE),
    sourceContentHash(0),
    sourceFaceCount(0),
    sourcePointCount(0),
    localToSource(SbMatrix::identity()),
    meshGeometry(FALSE),
    lodBacked(FALSE),
    sourceMeshRequestValid(FALSE),
    visible(TRUE),
    selected(FALSE),
    highlighted(FALSE)
{
    meshAssetBounds.makeEmpty();
    localBounds.makeEmpty();
}

BObolCompactLodPlanningSummary::BObolCompactLodPlanningSummary(void) :
    valid(FALSE),
    localToSource(SbMatrix::identity()),
    meshGeometry(FALSE),
    lodBacked(FALSE),
    sourceMeshRequestValid(FALSE),
    visible(TRUE),
    selected(FALSE),
    highlighted(FALSE)
{
    localBounds.makeEmpty();
}

BObolExternalLineSet::BObolExternalLineSet(void) :
    points(NULL),
    commands(NULL),
    precisePoints(NULL),
    count(0),
    sourceType("line-set"),
    geometryKind("line")
{
}

BObolExternalPointSet::BObolExternalPointSet(void) :
    points(NULL),
    precisePoints(NULL),
    count(0),
    sourceType("point-set"),
    geometryKind("point")
{
}

BObolExternalTriangleMesh::BObolExternalTriangleMesh(void) :
    points(NULL),
    pointCount(0),
    indices(NULL),
    indexCount(0),
    normals(NULL),
    normalCount(0),
    sourceType("indexed-face-set"),
    geometryKind("surface"),
    lodBacked(FALSE)
{
}

BObolExternalAnnotationSegment::BObolExternalAnnotationSegment(void) :
    kind(SEGMENT_NONE),
    lineStart(0),
    lineEnd(0),
    textRefPoint(0),
    text(NULL)
{
}

BObolExternalAnnotation::BObolExternalAnnotation(void) :
    basePoint(0.0f, 0.0f, 0.0f),
    linePoints(NULL),
    lineCommands(NULL),
    preciseLinePoints(NULL),
    linePointCount(0),
    annotationPoints(NULL),
    preciseAnnotationPoints(NULL),
    annotationPointCount(0),
    segments(NULL),
    segmentCount(0),
    sourceType("annotation"),
    geometryKind("annotation")
{
}

BObolDatabaseSourceDisplayPatch::BObolDatabaseSourceDisplayPatch(void) :
    visibleValid(FALSE),
    visible(TRUE),
    selectedValid(FALSE),
    selected(FALSE),
    highlightedValid(FALSE),
    highlighted(FALSE),
    lineStyleValid(FALSE),
    lineStyle(0),
    lineWidthValid(FALSE),
    lineWidth(0),
    transparencyValid(FALSE),
    transparency(0.0f),
    colorOverrideValid(FALSE),
    colorOverride(FALSE),
    colorValid(FALSE),
    color(1.0f, 1.0f, 1.0f),
    selectedColorValid(FALSE),
    selectedColor(1.0f, 1.0f, 1.0f),
    highlightedColorValid(FALSE),
    highlightedColor(1.0f, 1.0f, 0.0f),
    ghostedColorValid(FALSE),
    ghostedColor(0.55f, 0.55f, 0.55f)
{
}

BObolRealizedShapeSummary::BObolRealizedShapeSummary(void) :
    valid(FALSE),
    shapeKind(SHAPE_UNKNOWN),
    path(""),
    sourceName(""),
    sourceType(""),
    sourceId(0),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    ownerDrawMode(0),
    ownerSourceRevision(0),
    ownerInputsRevision(0),
    ownerViewRevision(0),
    ownerRealizedRevision(0),
    ownerRealizedSourceRevision(0),
    ownerRealizedInputsRevision(0),
    ownerRealizedViewRevision(0),
    ownerRealizationStatus(SoBRLDatabaseSource::UNREALIZED),
    ownerRealizationDiagnostic(""),
    ownerRealizationIdentity(""),
    ownerSourceStale(FALSE),
    ownerStaleReason(SoBRLDatabaseSource::STALE_NONE),
    displayName(""),
    geometryName(""),
    cacheIdentity(""),
    sourceIdentity(""),
    databaseIntent(FALSE),
    overlayIntent(FALSE),
    hudIntent(FALSE),
    localSource(FALSE),
    sharedSource(FALSE),
    nonDatabaseSource(FALSE),
    drawMode(0),
    recordRole(""),
    geometryKind(""),
    regionId(0),
    airCode(0),
    materialId(0),
    los(0),
    materialColorValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialShader(""),
    materialRevision(0),
    visible(FALSE),
    selectable(FALSE),
    selected(FALSE),
    highlighted(FALSE),
    ghosted(FALSE),
    hiddenLine(FALSE),
    editEmphasis(FALSE),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0f),
    editIntentId(""),
    editIntentRole(""),
    lodPolicy(0),
    colorOverride(FALSE),
    color(1.0f, 1.0f, 1.0f),
    lodAvailable(FALSE),
    lodActiveLevel(-1),
    lodFaceCount(0),
    lodPointCount(0),
    lodOriginalPointCount(0),
    lodNormalCount(0),
    lodHasSnappedPoints(FALSE),
    lodHasNormals(FALSE),
    lodBoundsMin(0.0f, 0.0f, 0.0f),
    lodBoundsMax(0.0f, 0.0f, 0.0f),
    pointCount(0),
    commandCount(0),
    segmentCount(0),
    pointPrimitiveCount(0),
    triangleCount(0),
    indexCount(0),
    boundsValid(FALSE),
    bounds()
{
    this->bounds.makeEmpty();
}

BObolCompactOccurrence::BObolCompactOccurrence(void) :
    geometry(),
    summary(),
    geometryTransform(SbMatrix::identity()),
    localTransform(SbMatrix::identity()),
    lodBacked(FALSE),
    sourceMeshRequestValid(FALSE),
    sourceMeshRequest(),
    occurrenceIndex(0),
    booleanOperation(SoBRLDatabaseSource::BOOLEAN_UNION)
{
}

BObolRealizedMaterialSummary::BObolRealizedMaterialSummary(void) :
    valid(FALSE),
    sourcePath(""),
    sourceName(""),
    sourceType(""),
    sourceId(0),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    ownerDrawMode(0),
    ownerSourceRevision(0),
    ownerInputsRevision(0),
    ownerViewRevision(0),
    ownerRealizedRevision(0),
    ownerRealizedSourceRevision(0),
    ownerRealizedInputsRevision(0),
    ownerRealizedViewRevision(0),
    ownerRealizationStatus(SoBRLDatabaseSource::UNREALIZED),
    ownerRealizationDiagnostic(""),
    ownerRealizationIdentity(""),
    ownerSourceStale(FALSE),
    ownerStaleReason(SoBRLDatabaseSource::STALE_NONE),
    materialName(""),
    parentName(""),
    materialSource(""),
    propertyCount(0)
{
}

BObolSceneTreeSummary::BObolSceneTreeSummary(void) :
    valid(FALSE),
    nodeKind(NODE_UNKNOWN),
    isGroup(FALSE),
    isShape(FALSE),
    isDatabaseSource(FALSE),
    isMaterialObject(FALSE),
    hasParent(FALSE),
    drawTreeDepth(0),
    childCount(0),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    path(""),
    sourceName(""),
    sourceType(""),
    sourceId(0),
    displayName(""),
    geometryName("")
{
}

BObolSceneDisplaySummary::BObolSceneDisplaySummary(void) :
    valid(FALSE),
    nodeKind(BObolSceneTreeSummary::NODE_UNKNOWN),
    isDatabaseSource(FALSE),
    hasDrawIntent(FALSE),
    intentPath(""),
    intentDrawMode(-1),
    visible(TRUE),
    selected(FALSE),
    highlighted(FALSE),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0),
    drawMode(0),
    materialValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialRevision(0),
    drawMatrixValid(FALSE),
    drawMatrix(SbMatrix::identity()),
    drawCenterValid(FALSE),
    drawCenter(0.0f, 0.0f, 0.0f),
    drawSizeValid(FALSE),
    drawSize(0.0f),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    path("")
{
}

BObolSceneMaterialSummary::BObolSceneMaterialSummary(void) :
    valid(FALSE),
    nodeKind(BObolSceneTreeSummary::NODE_UNKNOWN),
    materialValid(FALSE),
    materialRevision(0),
    materialColor(1.0f, 1.0f, 1.0f),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    path("")
{
}

BObolSceneBoundsSummary::BObolSceneBoundsSummary(void) :
    valid(FALSE),
    nodeKind(BObolSceneTreeSummary::NODE_UNKNOWN),
    boundsValid(FALSE),
    bounds(),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    path("")
{
    this->bounds.makeEmpty();
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

static const char *
database_source_skip_leading_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

static SbString
record_identity_with_revision(const char *identity, uint32_t revision)
{
    SbString ret = identity ? identity : "";
    char revisionString[64] = {0};
    snprintf(revisionString, sizeof(revisionString), "#%u", revision);
    ret += revisionString;
    return ret;
}

static int
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

static SbString
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
    uint64_t hash = 1469598103934665603ULL;
    const uintptr_t databaseIdentity = reinterpret_cast<uintptr_t>(
	source ? source->getDatabase() : NULL);
    const unsigned char *dbBytes = reinterpret_cast<const unsigned char *>(
	&databaseIdentity);
    for (size_t i = 0; i < sizeof(databaseIdentity); i++) {
	hash ^= dbBytes[i];
	hash *= 1099511628211ULL;
    }
    const SbString key = source_effective_instance_key(source);
    const unsigned char *bytes = reinterpret_cast<const unsigned char *>(
	key.getString());
    for (int i = 0; i < key.getLength(); i++) {
	hash ^= bytes[i];
	hash *= 1099511628211ULL;
    }
    return hash ? hash : 1;
}

static SbString
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
source_realization_identity(const SoBRLDatabaseSource *source)
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
		     source->sourceRevision.getValue(),
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
	if (!points.empty())
	    bobol_performance_counter_add(BOBOL_PERF_VLIST_POINTS,
		static_cast<uint64_t>(points.size()));
    }
    RT_FREE_VLIST(&rt_vlfree, &vhead);
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
	const std::vector<int32_t> &commands, Obol::PartGeometry &geometry)
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
	Obol::PartGeometry &geometry)
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

static int
cad_wire_part_geometry_from_bot(const struct rt_bot_internal *bot,
	Obol::PartGeometry &geometry)
{
    if (!bot || !bot->vertices || !bot->faces || !bot->num_vertices ||
	!bot->num_faces ||
	bot->num_faces > std::numeric_limits<size_t>::max() / 6u ||
	bot->num_faces > std::numeric_limits<uint32_t>::max() / 3u)
	return 0;
    RT_BOT_CK_MAGIC(bot);

    Obol::WireRep wire;
    wire.bounds.makeEmpty();
    wire.segmentPoints.reserve(bot->num_faces * 6u);
    wire.segmentIds.reserve(bot->num_faces * 3u);
    uint32_t segmentId = 0;
    for (size_t i = 0; i < bot->num_faces; i++) {
	const int *face = &bot->faces[i * 3];
	if (face[0] < 0 || face[1] < 0 || face[2] < 0 ||
	    static_cast<size_t>(face[0]) >= bot->num_vertices ||
	    static_cast<size_t>(face[1]) >= bot->num_vertices ||
	    static_cast<size_t>(face[2]) >= bot->num_vertices)
	    continue;
	const int edge[3][2] = {
	    {face[0], face[1]}, {face[1], face[2]}, {face[2], face[0]}
	};
	for (size_t j = 0; j < 3; j++) {
	    const fastf_t *a = &bot->vertices[edge[j][0] * 3];
	    const fastf_t *b = &bot->vertices[edge[j][1] * 3];
	    const SbVec3f start(static_cast<float>(a[X]),
		static_cast<float>(a[Y]), static_cast<float>(a[Z]));
	    const SbVec3f end(static_cast<float>(b[X]),
		static_cast<float>(b[Y]), static_cast<float>(b[Z]));
	    wire.segmentPoints.push_back(start);
	    wire.segmentPoints.push_back(end);
	    wire.segmentIds.push_back(segmentId++);
	    wire.bounds.extendBy(start);
	    wire.bounds.extendBy(end);
	}
    }
    if (wire.segmentPoints.empty())
	return 0;
    geometry.wire = std::move(wire);
    bobol_performance_counter_add(BOBOL_PERF_VLIST_POINTS,
	static_cast<uint64_t>(bot->num_faces) * 4u);
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

static SbBool
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

static void
remove_non_auxiliary_children(SoGroup *group)
{
    if (!group)
	return;

    for (int i = group->getNumChildren() - 1; i >= 0; i--) {
	SoNode *child = group->getChild(i);
	if (node_is_source_placement_transform(child) ||
	    node_is_auxiliary_vlist(child) ||
	    node_is_auxiliary_source(child))
	    continue;
	group->removeChild(i);
    }
}

template <typename ShapeT>
static void
unref_cached_geometry_map(std::map<std::string, ShapeT *> &cached)
{
    for (typename std::map<std::string, ShapeT *>::iterator it = cached.begin();
	 it != cached.end(); ++it) {
	if (it->second)
	    it->second->unref();
    }
    cached.clear();
}

template <typename ShapeT>
static void
store_cached_geometry_map(std::map<std::string, ShapeT *> &cached,
			  const std::string &key,
			  ShapeT *shape)
{
    if (!shape)
	return;

    /*
     * Realization seeding rebuilds the per-action cache from already realized
     * shared geometry.  Keep the internal shared-geometry key aligned with the
     * exact cache key used for storage, including view-LoD policy suffixes.
     * Instance shapes keep their user-facing geometry names separately.
     */
    shape->geometryName = key.c_str();
    shape->cacheIdentity =
	record_identity_with_revision(key.c_str(), shape->sourceId.getValue());

    typename std::map<std::string, ShapeT *>::iterator found =
	cached.find(key);
    if (found != cached.end()) {
	if (found->second == shape)
	    return;
	if (found->second)
	    found->second->unref();
    }

    shape->ref();
    cached[key] = shape;
}

BObolDatabaseSourceRealizationCache::BObolDatabaseSourceRealizationCache(void) :
    preserveCompactSourceOnFailure(false)
{
}

BObolDatabaseSourceRealizationCache::~BObolDatabaseSourceRealizationCache(void)
{
    this->clear();
}

void
BObolDatabaseSourceRealizationCache::clear(void)
{
    unref_cached_geometry_map(this->sharedWireGeometry);
    unref_cached_geometry_map(this->sharedMeshVListGeometry);
    unref_cached_geometry_map(this->sharedMeshGeometry);
    this->sharedWireBounds.clear();
    this->sharedWireCadGeometry.clear();
    this->sharedMeshVListCadGeometry.clear();
    this->sharedMeshCadGeometry.clear();
}

template <typename ValueT>
static void
erase_realization_cache_values(std::map<std::string, ValueT> &values,
	const std::string &name)
{
    values.erase(name);
    const std::string prefix = name + ":";
    for (auto it = values.lower_bound(prefix); it != values.end() &&
	it->first.compare(0, prefix.size(), prefix) == 0;)
	it = values.erase(it);
}

template <typename ShapeT>
static void
erase_realization_cache_shapes(std::map<std::string, ShapeT *> &values,
	const std::string &name)
{

    auto exact = values.find(name);
    if (exact != values.end()) {
	if (exact->second)
	    exact->second->unref();
	values.erase(exact);
    }
    const std::string prefix = name + ":";
    for (auto it = values.lower_bound(prefix); it != values.end() &&
	it->first.compare(0, prefix.size(), prefix) == 0;) {
	if (it->second)
	    it->second->unref();
	it = values.erase(it);
    }
}

void
BObolDatabaseSourceRealizationCache::eraseObject(const std::string &path)
{
    const size_t slash = path.find_last_of('/');
    const std::string name = slash == std::string::npos ? path :
	path.substr(slash + 1);
    if (name.empty())
	return;
    erase_realization_cache_shapes(this->sharedWireGeometry, name);
    erase_realization_cache_values(this->sharedWireBounds, name);
    erase_realization_cache_shapes(this->sharedMeshVListGeometry, name);
    erase_realization_cache_shapes(this->sharedMeshGeometry, name);
    erase_realization_cache_values(this->sharedWireCadGeometry, name);
    erase_realization_cache_values(this->sharedMeshVListCadGeometry, name);
    erase_realization_cache_values(this->sharedMeshCadGeometry, name);
}

static std::string
realization_cache_object_name(const std::string &path)
{
    const size_t slash = path.find_last_of('/');
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

template <typename ValueT>
static void
rename_realization_cache_values(std::map<std::string, ValueT> &values,
	const std::string &oldName, const std::string &newName)
{
    std::vector<std::string> keys;
    if (values.find(oldName) != values.end())
	keys.push_back(oldName);
    const std::string prefix = oldName + ":";
    for (auto it = values.lower_bound(prefix); it != values.end() &&
	it->first.compare(0, prefix.size(), prefix) == 0; ++it)
	keys.push_back(it->first);
    for (const std::string &key : keys) {
	auto it = values.find(key);
	if (it == values.end())
	    continue;
	const std::string nextKey = newName + key.substr(oldName.size());
	ValueT value = std::move(it->second);
	values.erase(it);
	values.erase(nextKey);
	values.emplace(nextKey, std::move(value));
    }
}

template <typename ShapeT>
static void
rename_realization_cache_shapes(std::map<std::string, ShapeT *> &values,
	const std::string &oldName, const std::string &newName)
{
    std::vector<std::string> keys;
    if (values.find(oldName) != values.end())
	keys.push_back(oldName);
    const std::string prefix = oldName + ":";
    for (auto it = values.lower_bound(prefix); it != values.end() &&
	it->first.compare(0, prefix.size(), prefix) == 0; ++it)
	keys.push_back(it->first);
    for (const std::string &key : keys) {
	auto it = values.find(key);
	if (it == values.end())
	    continue;
	const std::string nextKey = newName + key.substr(oldName.size());
	ShapeT *shape = it->second;
	values.erase(it);
	auto destination = values.find(nextKey);
	if (destination != values.end()) {
	    if (destination->second)
		destination->second->unref();
	    values.erase(destination);
	}
	values.emplace(nextKey, shape);
    }
}

void
BObolDatabaseSourceRealizationCache::renameObject(
    const std::string &oldPath, const std::string &newPath)
{
    const std::string oldName = realization_cache_object_name(oldPath);
    const std::string newName = realization_cache_object_name(newPath);
    if (oldName.empty() || newName.empty() || oldName == newName)
	return;
    rename_realization_cache_shapes(this->sharedWireGeometry, oldName,
	newName);
    rename_realization_cache_values(this->sharedWireBounds, oldName, newName);
    rename_realization_cache_shapes(this->sharedMeshVListGeometry, oldName,
	newName);
    rename_realization_cache_shapes(this->sharedMeshGeometry, oldName,
	newName);
    rename_realization_cache_values(this->sharedWireCadGeometry, oldName,
	newName);
    rename_realization_cache_values(this->sharedMeshVListCadGeometry, oldName,
	newName);
    rename_realization_cache_values(this->sharedMeshCadGeometry, oldName,
	newName);
}

template <typename ValueT>
static void
erase_realization_cache_variant_values(
    std::map<std::string, ValueT> &values, const char *marker)
{
    for (auto it = values.begin(); it != values.end();) {
	if (it->first.find(marker) != std::string::npos)
	    it = values.erase(it);
	else
	    ++it;
    }
}

template <typename ShapeT>
static void
erase_realization_cache_variant_shapes(
    std::map<std::string, ShapeT *> &values, const char *marker)
{
    for (auto it = values.begin(); it != values.end();) {
	if (it->first.find(marker) == std::string::npos) {
	    ++it;
	    continue;
	}
	if (it->second)
	    it->second->unref();
	it = values.erase(it);
    }
}

void
BObolDatabaseSourceRealizationCache::eraseViewVariants(void)
{
    static const char marker[] = ":view-lod:";
    erase_realization_cache_variant_shapes(this->sharedWireGeometry, marker);
    erase_realization_cache_variant_values(this->sharedWireBounds, marker);
    erase_realization_cache_variant_shapes(this->sharedMeshVListGeometry,
	marker);
    erase_realization_cache_variant_shapes(this->sharedMeshGeometry, marker);
    erase_realization_cache_variant_values(this->sharedWireCadGeometry,
	marker);
    erase_realization_cache_variant_values(this->sharedMeshVListCadGeometry,
	marker);
    erase_realization_cache_variant_values(this->sharedMeshCadGeometry,
	marker);
}

void
BObolDatabaseSourceRealizationCache::storeWireGeometry(
    const std::string &key,
    SoBRLVListShape *shape)
{
    store_cached_geometry_map(this->sharedWireGeometry, key, shape);
}

void
BObolDatabaseSourceRealizationCache::storeWireBounds(
    const std::string &key,
    const SbBox3f &bounds)
{
    if (key.empty() || bounds.isEmpty())
	return;

    this->sharedWireBounds[key] = bounds;
}

void
BObolDatabaseSourceRealizationCache::storeMeshVListGeometry(
    const std::string &key,
    SoBRLVListShape *shape)
{
    store_cached_geometry_map(this->sharedMeshVListGeometry, key, shape);
}

void
BObolDatabaseSourceRealizationCache::storeMeshGeometry(
    const std::string &key,
    SoBRLMeshShape *shape)
{
    store_cached_geometry_map(this->sharedMeshGeometry, key, shape);
}

static std::shared_ptr<const Obol::PartGeometry>
store_cached_part_geometry(
    std::map<std::string, BObolCachedPartGeometry> &cache,
    const std::string &key, Obol::PartGeometry &&geometry,
    const char *sourceType, const char *geometryKind, const SbBox3f *bounds,
    bool lodBacked, const BObolSourceMeshRequest *sourceMeshRequest)
{
    if (key.empty())
	return std::shared_ptr<const Obol::PartGeometry>();

    BObolCachedPartGeometry &stored = cache[key];
    stored.geometry =
	std::make_shared<const Obol::PartGeometry>(std::move(geometry));
    stored.sourceType = sourceType ? sourceType : "";
    stored.geometryKind = geometryKind ? geometryKind : "";
    if (bounds)
	stored.bounds = *bounds;
    else
	stored.bounds = compact_part_geometry_bounds(stored.geometry);
    stored.geometryTransform = SbMatrix::identity();
    stored.lodBacked = lodBacked;
    stored.sourceMeshRequestValid = sourceMeshRequest != NULL;
    if (sourceMeshRequest)
	stored.sourceMeshRequest = *sourceMeshRequest;
    return stored.geometry;
}

static std::shared_ptr<const Obol::PartGeometry>
store_cached_part_geometry_reference(
    std::map<std::string, BObolCachedPartGeometry> &cache,
    const std::string &key, const std::shared_ptr<const Obol::PartGeometry> &geometry,
    const SbMatrix &geometryTransform, const char *sourceType,
    const char *geometryKind, const SbBox3f *bounds, bool lodBacked,
    const BObolSourceMeshRequest *sourceMeshRequest)
{
    if (key.empty() || !geometry)
	return std::shared_ptr<const Obol::PartGeometry>();

    BObolCachedPartGeometry &stored = cache[key];
    stored.geometry = geometry;
    stored.sourceType = sourceType ? sourceType : "";
    stored.geometryKind = geometryKind ? geometryKind : "";
    if (bounds)
	stored.bounds = *bounds;
    else
	stored.bounds = database_source_transform_bounds(
	    compact_part_geometry_bounds(geometry), geometryTransform);
    stored.geometryTransform = geometryTransform;
    stored.lodBacked = lodBacked;
    stored.sourceMeshRequestValid = sourceMeshRequest != NULL;
    if (sourceMeshRequest)
	stored.sourceMeshRequest = *sourceMeshRequest;
    return stored.geometry;
}

std::shared_ptr<const Obol::PartGeometry>
BObolDatabaseSourceRealizationCache::storeWireCadGeometry(
    const std::string &key, Obol::PartGeometry &&geometry,
    const char *sourceType, const char *geometryKind, const SbBox3f *bounds,
    bool lodBacked, const BObolSourceMeshRequest *sourceMeshRequest)
{
    return store_cached_part_geometry(this->sharedWireCadGeometry, key,
	std::move(geometry), sourceType, geometryKind, bounds, lodBacked,
	sourceMeshRequest);
}

std::shared_ptr<const Obol::PartGeometry>
BObolDatabaseSourceRealizationCache::storeMeshVListCadGeometry(
    const std::string &key, Obol::PartGeometry &&geometry,
    const char *sourceType, const char *geometryKind, const SbBox3f *bounds,
    bool lodBacked, const BObolSourceMeshRequest *sourceMeshRequest)
{
    return store_cached_part_geometry(this->sharedMeshVListCadGeometry, key,
	std::move(geometry), sourceType, geometryKind, bounds, lodBacked,
	sourceMeshRequest);
}

std::shared_ptr<const Obol::PartGeometry>
BObolDatabaseSourceRealizationCache::storeMeshCadGeometry(
    const std::string &key, Obol::PartGeometry &&geometry,
    const char *sourceType, const char *geometryKind, const SbBox3f *bounds,
    bool lodBacked, const BObolSourceMeshRequest *sourceMeshRequest)
{
    return store_cached_part_geometry(this->sharedMeshCadGeometry, key,
	std::move(geometry), sourceType, geometryKind, bounds, lodBacked,
	sourceMeshRequest);
}

std::shared_ptr<const Obol::PartGeometry>
BObolDatabaseSourceRealizationCache::storeMeshCadGeometryReference(
    const std::string &key,
    const std::shared_ptr<const Obol::PartGeometry> &geometry,
    const SbMatrix &geometryTransform, const char *sourceType,
    const char *geometryKind, const SbBox3f *bounds, bool lodBacked,
    const BObolSourceMeshRequest *sourceMeshRequest)
{
    return store_cached_part_geometry_reference(this->sharedMeshCadGeometry,
	key, geometry, geometryTransform, sourceType, geometryKind, bounds,
	lodBacked, sourceMeshRequest);
}

static const BObolCachedPartGeometry *
find_cached_part_geometry(
    const std::map<std::string, BObolCachedPartGeometry> &cache,
    const std::string &key)
{
    auto found = cache.find(key);
    return found == cache.end() || !found->second.geometry ? NULL :
	&found->second;
}

const BObolCachedPartGeometry *
BObolDatabaseSourceRealizationCache::findWireCadGeometry(
    const std::string &key) const
{
    return find_cached_part_geometry(this->sharedWireCadGeometry, key);
}

const BObolCachedPartGeometry *
BObolDatabaseSourceRealizationCache::findMeshVListCadGeometry(
    const std::string &key) const
{
    return find_cached_part_geometry(this->sharedMeshVListCadGeometry, key);
}

const BObolCachedPartGeometry *
BObolDatabaseSourceRealizationCache::findMeshCadGeometry(
    const std::string &key) const
{
    return find_cached_part_geometry(this->sharedMeshCadGeometry, key);
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

    std::map<std::string, BObolCachedPartGeometry>::iterator found =
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
		cache->storeMeshVListGeometry(key, geometry);
	    else
		cache->storeWireGeometry(key, geometry);
	}
	return;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	SoBRLMeshShape *geometry = shape->getSharedGeometrySource();
	std::string key;
	if (meshRealization && cacheable_shared_geometry_key(shape,
		geometry, key))
	    cache->storeMeshGeometry(key, geometry);
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
	compact_index(NULL),
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
    BObolCompactInstanceIndex *compact_index;
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

static int cad_vlist_part_geometry(const SoBRLVListShape *shape,
	Obol::PartGeometry &geometry);
static int cad_mesh_part_geometry(const SoBRLMeshShape *shape,
	Obol::PartGeometry &geometry);
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
	int &unsupported);
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
	    Obol::CadIdBuilder::hash128(occurrenceKey);
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

static SbMatrix
mat_to_sbmatrix(const mat_t mat)
{
    return SbMatrix(
	       static_cast<float>(mat[0]), static_cast<float>(mat[4]),
	       static_cast<float>(mat[8]), static_cast<float>(mat[12]),
	       static_cast<float>(mat[1]), static_cast<float>(mat[5]),
	       static_cast<float>(mat[9]), static_cast<float>(mat[13]),
	       static_cast<float>(mat[2]), static_cast<float>(mat[6]),
	       static_cast<float>(mat[10]), static_cast<float>(mat[14]),
	       static_cast<float>(mat[3]), static_cast<float>(mat[7]),
	       static_cast<float>(mat[11]), static_cast<float>(mat[15]));
}

static SoSeparator *
realize_instance_leaf_separator(const struct db_tree_state *tsp)
{
    SoSeparator *leaf = new SoSeparator;
    if (!tsp || bn_mat_is_identity(tsp->ts_mat))
	return leaf;

    SoMatrixTransform *transform = new SoMatrixTransform;
    transform->matrix = mat_to_sbmatrix(tsp->ts_mat);
    leaf->addChild(transform);
    return leaf;
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
			 const SoBRLDatabaseSource *source)
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
    summary.materialColorValid = source ?
	source->materialColorValid.getValue() : FALSE;
    summary.materialColor = source ? source->materialColor.getValue() :
	SbColor(1.0f, 1.0f, 1.0f);
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
	const struct db_tree_state *tsp, const char *path,
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
    if (source && source->materialColorValid.getValue() &&
	(source->materialPolicy.getValue() !=
	 SoBRLDatabaseSource::MATERIAL_DATABASE ||
	 !summary.materialColorValid)) {
	summary.materialColorValid = TRUE;
	summary.materialColor = source->materialColor.getValue();
    }
    /* db_walk_tree has already resolved combination inheritance into tsp.
     * Re-walking the full database path for every leaf dominated large-model
     * realization.  Only consult db_full_path_color when the tree state has no
     * color (notably region-id color-table fallback). */
    if (source && !source->colorOverride.getValue() &&
	!summary.materialColorValid) {
	BObolMaterialPathState resolved;
	SbColor databaseColor;
	const bool haveColor = materialSweep ?
	    materialSweep->resolve(path, resolved) :
	    bobol_database_source_path_material_color(source->getDatabase(),
		path, databaseColor);
	if (haveColor) {
	    summary.materialColorValid = TRUE;
	    summary.materialColor = materialSweep ? resolved.color :
		databaseColor;
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
    summary.lodActiveLevel = request.lodActiveLevel;
    summary.lodFaceCount = request.lodFaceCount;
    summary.lodPointCount = request.lodPointCount;
    summary.lodOriginalPointCount = request.lodOriginalPointCount;
    summary.lodNormalCount = request.lodNormalCount;
    summary.lodHasSnappedPoints = request.lodHasSnappedPoints ? TRUE : FALSE;
    summary.lodHasNormals = request.lodHasNormals ? TRUE : FALSE;
    summary.lodBoundsMin = request.lodBoundsMin;
    summary.lodBoundsMax = request.lodBoundsMax;
}

static void
retarget_realized_node_source(SoNode *node,
			      const SoBRLDatabaseSource *source,
			      const char *oldSourcePath,
			      const char *newSourcePath,
			      uint32_t revision);

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

    std::string nextDisplayName;
    if (source->displayName.getValue().getLength() > 0) {
	nextDisplayName = source->displayName.getValue().getString();
    } else if (shape->sourceName.getValue().getLength() > 0) {
	nextDisplayName = shape->sourceName.getValue().getString();
    } else {
	nextDisplayName = shape->sourcePath.getValue().getString();
    }
    database_source_assign_string(shape->displayName,
				  nextDisplayName.c_str());
}

static void
sync_realized_shape_owner_state_in_node(SoNode *node,
					const SoBRLDatabaseSource *source)
{
    if (!node || !source)
	return;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	sync_shape_owner_state(shape, source);
	sync_shape_placement_state(shape, source);
	sync_shape_display_state(shape, source);
	sync_shape_display_name(shape, source);
	return;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	sync_shape_owner_state(shape, source);
	sync_shape_placement_state(shape, source);
	sync_shape_display_state(shape, source);
	sync_shape_display_name(shape, source);
	return;
    }
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    sync_realized_shape_owner_state_in_node(group->getChild(i), source);
    }
}

static void
sync_realized_shape_placement_state_in_node(SoNode *node,
				    const SoBRLDatabaseSource *source)
{
    if (!node || !source)
	return;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	sync_shape_placement_state(static_cast<SoBRLVListShape *>(node),
		source);
	return;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	sync_shape_placement_state(static_cast<SoBRLMeshShape *>(node),
		source);
	return;
    }
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    sync_realized_shape_placement_state_in_node(group->getChild(i),
		    source);
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

    point_t annotationBase = VINIT_ZERO;
    struct rt_annot_internal *annotation = NULL;
    const char *typeLabel = primitive_type_label(intern);
    const bool addAnnotationBase = primitive_is_annotation(intern->idb_type,
				   typeLabel);
    if (addAnnotationBase) {
	struct rt_annot_internal *annot =
		static_cast<struct rt_annot_internal *>(intern->idb_ptr);
	RT_ANNOT_CK_MAGIC(annot);
	annotation = rt_copy_annot(annot);
	VMOVE(annotationBase, annot->V);
    }

    struct bg_tess_tol ttol = plotTtol ? *plotTtol : source_tess_tol(source);
    struct bn_tol tol = BN_TOL_INIT_TOL;
    if (plotTol)
	tol = *plotTol;
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    const int ret = plot_internal_to_vlist_geometry(points, commands,
	intern, &ttol, &tol);
    if (ret < 0) {
	free_annotation_record_copy(annotation);
	return NULL;
    }

    if (addAnnotationBase) {
	for (size_t i = 0; i < points.size(); i++) {
	    points[i].setValue(points[i][0] +
			       static_cast<float>(annotationBase[X]),
			       points[i][1] + static_cast<float>(annotationBase[Y]),
			       points[i][2] + static_cast<float>(annotationBase[Z]));
	}
    }
    if (points.empty() || points.size() != commands.size()) {
	free_annotation_record_copy(annotation);
	return NULL;
    }

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->setLineSet(points.data(), commands.data(),
		      static_cast<int>(points.size()));
    assign_annotation_record(shape, annotation);
    free_annotation_record_copy(annotation);
    return shape;
}

static int
cad_wire_part_geometry_from_plot_internal(struct rt_db_internal *intern,
	const SoBRLDatabaseSource *source, Obol::PartGeometry &geometry)
{
    if (!internal_payload_magic_valid(intern) || !intern->idb_meth ||
	!intern->idb_meth->ft_plot)
	return 0;

    point_t annotationBase = VINIT_ZERO;
    const char *typeLabel = primitive_type_label(intern);
    const bool addAnnotationBase = primitive_is_annotation(intern->idb_type,
	typeLabel);
    if (addAnnotationBase) {
	struct rt_annot_internal *annotation =
	    static_cast<struct rt_annot_internal *>(intern->idb_ptr);
	RT_ANNOT_CK_MAGIC(annotation);
	VMOVE(annotationBase, annotation->V);
    }

    struct bg_tess_tol ttol = source_tess_tol(source);
    struct bn_tol tol = BN_TOL_INIT_TOL;
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    if (plot_internal_to_vlist_geometry(points, commands, intern,
	    &ttol, &tol) < 0)
	return 0;
    if (addAnnotationBase) {
	for (SbVec3f &point : points) {
	    point.setValue(point[0] + static_cast<float>(annotationBase[X]),
		point[1] + static_cast<float>(annotationBase[Y]),
		point[2] + static_cast<float>(annotationBase[Z]));
	}
    }
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
	double tolerance)
{
    std::vector<size_t> result;
    if (points.size() < 3 || !(tolerance > 0.0)) {
	result.resize(points.size());
	std::iota(result.begin(), result.end(), 0);
	return result;
    }
    std::vector<unsigned char> keep(points.size(), 0);
    keep.front() = 1;
    keep.back() = 1;
    std::vector<std::pair<size_t, size_t>> work;
    work.emplace_back(0, points.size() - 1);
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
	    keep[split] = 1;
	    if (split > range.first + 1)
		work.emplace_back(range.first, split);
	    if (range.second > split + 1)
		work.emplace_back(split, range.second);
	}
    }
    result.reserve(points.size());
    for (size_t i = 0; i < keep.size(); ++i)
	if (keep[i])
	    result.push_back(i);
    return result;
}

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
    Obol::PartGeometry &geometry)
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

    std::vector<std::vector<SbVec3f>> curves;
    for (size_t i = 0; i < realization.line_count; ++i) {
	const int command = realization.line_commands ?
	    realization.line_commands[i] : RT_PRIMITIVE_LINE_DRAW;
	const SbVec3f point(
	    static_cast<float>(realization.line_points[i][X]),
	    static_cast<float>(realization.line_points[i][Y]),
	    static_cast<float>(realization.line_points[i][Z]));
	if (command == RT_PRIMITIVE_LINE_MOVE || curves.empty())
	    curves.emplace_back();
	if (curves.back().empty() || curves.back().back() != point)
	    curves.back().push_back(point);
    }
    primitive_realization_line_set_free(&realization);
    curves.erase(std::remove_if(curves.begin(), curves.end(),
	[](const std::vector<SbVec3f> &curve) { return curve.size() < 2; }),
	curves.end());
    if (curves.empty())
	return 0;

    Obol::WireRep wire;
    wire.bounds.makeEmpty();
    for (const std::vector<SbVec3f> &curve : curves)
	for (const SbVec3f &point : curve)
	    wire.bounds.extendBy(point);
    if (wire.bounds.isEmpty())
	return 0;

    std::vector<SbVec3f> previousPoints;
    std::vector<uint32_t> previousIds;
    uint32_t previousFirst = 0;
    uint32_t previousCount = 0;
    for (uint8_t level = 0; level < 16; ++level) {
	std::vector<SbVec3f> levelPoints;
	std::vector<uint32_t> levelIds;
	for (size_t curveIndex = 0; curveIndex < curves.size(); ++curveIndex) {
	    const std::vector<SbVec3f> &curve = curves[curveIndex];
	    SbBox3f curveBounds;
	    curveBounds.makeEmpty();
	    for (const SbVec3f &point : curve)
		curveBounds.extendBy(point);
	    const SbVec3f extent = curveBounds.getMax() - curveBounds.getMin();
	    const double diagonal = sqrt(static_cast<double>(extent.sqrLength()));
	    const double tolerance = level == 15 ? 0.0 :
		ldexp(diagonal, -static_cast<int>(level) - 2);
	    const std::vector<size_t> selected =
		cad_wire_simplify_polyline(curve, tolerance);
	    for (size_t i = 1; i < selected.size(); ++i) {
		levelPoints.push_back(curve[selected[i - 1]]);
		levelPoints.push_back(curve[selected[i]]);
		levelIds.push_back(static_cast<uint32_t>(
		    std::min<size_t>(curveIndex, UINT32_MAX)));
	    }
	}
	if (levelPoints.size() / 2 > UINT32_MAX ||
	    wire.segmentCount() > UINT32_MAX - levelPoints.size() / 2)
	    return 0;
	if (level > 0 && cad_wire_level_equal(levelPoints, levelIds,
		previousPoints, previousIds)) {
	    wire.progressiveSegmentFirst[level] = previousFirst;
	    wire.progressiveSegmentCount[level] = previousCount;
	    continue;
	}
	const uint32_t first = static_cast<uint32_t>(wire.segmentCount());
	const uint32_t count = static_cast<uint32_t>(levelPoints.size() / 2);
	wire.segmentPoints.insert(wire.segmentPoints.end(),
	    levelPoints.begin(), levelPoints.end());
	wire.segmentIds.insert(wire.segmentIds.end(),
	    levelIds.begin(), levelIds.end());
	wire.progressiveSegmentFirst[level] = first;
	wire.progressiveSegmentCount[level] = count;
	previousPoints = std::move(levelPoints);
	previousIds = std::move(levelIds);
	previousFirst = first;
	previousCount = count;
    }
    wire.progressiveMinimumLevel = 0;
    wire.progressiveResidentLevel = 15;
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
	const SbBox3f &localBounds, Obol::PartGeometry &geometry)
{
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
    return cad_wire_part_geometry_from_line_set(points, commands, geometry);
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

	    Obol::PartGeometry generated;
	    int generatedGeometry = localIntern->idb_type == ID_BOT ?
		cad_wire_part_geometry_from_bot(
		    static_cast<const struct rt_bot_internal *>(localIntern->idb_ptr),
		    generated) : cad_wire_part_geometry_from_lod_realization_internal(
			localIntern, data->source, localBounds, generated);
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
		localBounds.isEmpty() ? NULL : &localBounds);
	}

	char *path = db_path_to_string(pathp);
	compact_occurrence_build input;
	input.occurrence.geometry = cadGeometry;
	input.occurrence.localTransform = mat_to_sbmatrix(tsp->ts_mat);
	input.occurrence.summary = compact_occurrence_tree_summary(
	    data->source, tsp, path, dp->d_namep,
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
	if (!bounds.isEmpty()) {
	    data->compact_bounds.extendBy(bounds);
	    data->compact_bounds_valid = TRUE;
	}
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
    std::map<std::string, SoBRLVListShape *>::iterator found =
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
	    Obol::PartGeometry generated;
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
		data->source, tsp, path, dp->d_namep,
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
	    if (!bounds.isEmpty()) {
		data->compact_bounds.extendBy(bounds);
		data->compact_bounds_valid = TRUE;
	    }
	}
    } else {
	BObolPerformanceTimer timer(BOBOL_PERF_REALIZED_INSTANCE_NODE_US);
	SoSeparator *leaf = realize_instance_leaf_separator(tsp);
	SoBRLVListShape *shape = new SoBRLVListShape;
	assign_realized_identity(shape, tsp, path, dp->d_namep, typeLabel,
	    data->revision, data->source);
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

static std::string
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

static void
set_source_bounds_from_local_box(SoBRLDatabaseSource *source,
				 const SbBox3f &bounds)
{
    if (!source)
	return;

    if (bounds.isEmpty()) {
	source->clearSourceBounds();
	return;
    }

    (void)source->setSourceBoundsState(TRUE, bounds.getMin(),
				       bounds.getMax());
}

static int
source_has_auxiliary_children(const SoBRLDatabaseSource *source);
static int
cad_vlist_part_geometry(const SoBRLVListShape *shape,
			Obol::PartGeometry &geometry);

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
    std::map<std::string, SoBRLVListShape *>::iterator found =
	cache->sharedWireGeometry.find(cacheKey);
    if (found != cache->sharedWireGeometry.end())
	sharedShape = found->second;
    bobol_performance_counter_add(sharedShape ? BOBOL_PERF_WIRE_CACHE_HITS :
	BOBOL_PERF_WIRE_CACHE_MISSES, 1);

    std::map<std::string, SbBox3f>::const_iterator boundsFound =
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
	    SoSeparator *leaf = new SoSeparator;
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
	SoSeparator *leaf = new SoSeparator;
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
    uint32_t revision)
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
	    Obol::PartGeometry generated;
	    if (cad_wire_part_geometry_from_aabb(lodProxyBounds, generated)) {
		geometryKind = "aabb";
		localBounds = lodProxyBounds;
		localBoundsValid = TRUE;
		cadGeometry = cache->storeWireCadGeometry(cacheKey,
		    std::move(generated), typeLabel, geometryKind, &localBounds,
		    true);
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
	    SoSeparator *leaf = new SoSeparator;
	    assign_material_identity(materialObject, fullPath.c_str(),
		dp->d_namep, typeLabel, revision);
	    leaf->addChild(materialObject);
	    database_source_add_realized_child(source, leaf);
	    source->clearSourceBounds();
	    return 1;
	}

	if (!cadGeometry) {
	    Obol::PartGeometry generated;
	    int generatedGeometry = 0;
	    if (validInternal.local.idb_type == ID_BOT)
		generatedGeometry = cad_wire_part_geometry_from_bot(
		    static_cast<const struct rt_bot_internal *>(
			validInternal.local.idb_ptr), generated);
	    else
		generatedGeometry =
		    cad_wire_part_geometry_from_lod_realization_internal(
			&validInternal.local, source, localBounds, generated);
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
		    localBoundsValid ? &localBounds : NULL);
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
    occurrence.summary = compact_occurrence_summary(source,
	fullPath.c_str(), dp->d_namep,
	geometryKind && BU_STR_EQUAL(geometryKind, "annotation") ?
	"annotation" : typeLabel,
	geometryKind, revision, BObolRealizedShapeSummary::SHAPE_VLIST);
    occurrence.occurrenceIndex = source->occurrenceIndex.getValue();
    occurrence.booleanOperation = source->booleanOperation.getValue();
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
	Obol::PartGeometry &geometry)
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

static void
sanitize_triangle_normals(std::vector<SbVec3f> &normals,
			  const std::vector<SbVec3f> &points,
			  const std::vector<int32_t> &triangles);

static int
canonicalize_corner_normal_mesh(Obol::TriMesh &mesh,
	const std::vector<SbVec3f> &cornerNormals);

/* Missing normals carry no author intent to smooth.  Preserve BRL-CAD's
 * default shaded-BoT semantics by using the geometric normal of each triangle
 * at all three corners.  The separate per-view smooth policy may explicitly
 * synthesize presentation normals later. */
template <typename Index>
static void
generate_flat_triangle_normals(std::vector<SbVec3f> &normals,
			       const std::vector<SbVec3f> &points,
			       const std::vector<Index> &triangles)
{
    normals.clear();
    if (points.empty() || triangles.empty() || triangles.size() % 3)
	return;

    const size_t faceCount = triangles.size() / 3;
    normals.reserve(triangles.size());
    for (size_t faceIndex = 0; faceIndex < faceCount; faceIndex++) {
	const size_t indexBase = faceIndex * 3;
	const uint64_t ia = static_cast<uint64_t>(triangles[indexBase]);
	const uint64_t ib = static_cast<uint64_t>(triangles[indexBase + 1]);
	const uint64_t ic = static_cast<uint64_t>(triangles[indexBase + 2]);
	if (ia >= points.size() || ib >= points.size() || ic >= points.size()) {
	    normals.clear();
	    return;
	}

	SbVec3f normal =
	    (points[static_cast<size_t>(ib)] - points[static_cast<size_t>(ia)]).cross(
		points[static_cast<size_t>(ic)] - points[static_cast<size_t>(ia)]);
	if (normal.length() > 0.0f)
	    normal.normalize();
	else
	    normal = SbVec3f(0.0f, 0.0f, 1.0f);
	normals.push_back(normal);
	normals.push_back(normal);
	normals.push_back(normal);
    }
}

template <typename Index>
static void
generate_smooth_triangle_normals(std::vector<SbVec3f> &normals,
				 const std::vector<SbVec3f> &points,
				 const std::vector<Index> &triangles,
				 float creaseAngleDegrees)
{
    normals.clear();
    if (points.empty() || triangles.empty() || triangles.size() % 3)
	return;

    const size_t faceCount = triangles.size() / 3;
    std::vector<SbVec3f> faceNormals(faceCount);
    std::vector<std::vector<size_t>> vertexFaces(points.size());
    for (size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex) {
	const size_t base = faceIndex * 3;
	const uint64_t ia = static_cast<uint64_t>(triangles[base]);
	const uint64_t ib = static_cast<uint64_t>(triangles[base + 1]);
	const uint64_t ic = static_cast<uint64_t>(triangles[base + 2]);
	if (ia >= points.size() || ib >= points.size() || ic >= points.size()) {
	    normals.clear();
	    return;
	}
	SbVec3f normal =
	    (points[static_cast<size_t>(ib)] -
	     points[static_cast<size_t>(ia)]).cross(
		points[static_cast<size_t>(ic)] -
		points[static_cast<size_t>(ia)]);
	if (normal.sqrLength() > 0.0f)
	    normal.normalize();
	else
	    normal.setValue(0.0f, 0.0f, 1.0f);
	faceNormals[faceIndex] = normal;
	vertexFaces[static_cast<size_t>(ia)].push_back(faceIndex);
	vertexFaces[static_cast<size_t>(ib)].push_back(faceIndex);
	vertexFaces[static_cast<size_t>(ic)].push_back(faceIndex);
    }

    const float clampedAngle =
	std::max(0.0f, std::min(180.0f, creaseAngleDegrees));
    const float creaseCosine =
	static_cast<float>(std::cos(clampedAngle * M_PI / 180.0));
    normals.reserve(triangles.size());
    for (size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex) {
	const SbVec3f &faceNormal = faceNormals[faceIndex];
	for (size_t corner = 0; corner < 3; ++corner) {
	    const size_t vertexIndex = static_cast<size_t>(
		triangles[faceIndex * 3 + corner]);
	    SbVec3f smoothNormal(0.0f, 0.0f, 0.0f);
	    for (const size_t adjacentFace : vertexFaces[vertexIndex]) {
		const SbVec3f &candidate = faceNormals[adjacentFace];
		if (candidate.dot(faceNormal) >= creaseCosine)
		    smoothNormal += candidate;
	    }
	    if (smoothNormal.sqrLength() > 0.0f)
		smoothNormal.normalize();
	    else
		smoothNormal = faceNormal;
	    normals.push_back(smoothNormal);
	}
    }
}

static bool
bot_authored_triangle_normals(std::vector<SbVec3f> &normals,
			     const struct rt_bot_internal *bot)
{
    normals.clear();
    if (!bot ||
	(bot->bot_flags &
	    (RT_BOT_HAS_SURFACE_NORMALS | RT_BOT_USE_NORMALS)) !=
	    (RT_BOT_HAS_SURFACE_NORMALS | RT_BOT_USE_NORMALS) ||
	!bot->normals || !bot->face_normals ||
	bot->num_face_normals < bot->num_faces)
	return false;

    normals.reserve(bot->num_faces * 3);
    for (size_t faceIndex = 0; faceIndex < bot->num_faces; faceIndex++) {
	for (size_t corner = 0; corner < 3; corner++) {
	    const size_t sourceCorner =
		(bot->orientation == RT_BOT_CW && corner > 0) ? 3 - corner : corner;
	    const int normalIndex =
		bot->face_normals[faceIndex * 3 + sourceCorner];
	    if (normalIndex < 0 ||
		static_cast<size_t>(normalIndex) >= bot->num_normals) {
		normals.clear();
		return false;
	    }
	    const fastf_t *normal = &bot->normals[static_cast<size_t>(normalIndex) * 3];
	    SbVec3f converted(static_cast<float>(normal[X]),
		static_cast<float>(normal[Y]), static_cast<float>(normal[Z]));
	    if (converted.length() <= 0.0f) {
		normals.clear();
		return false;
	    }
	    converted.normalize();
	    normals.push_back(converted);
	}
    }
    return true;
}

template <typename Index>
static void
bot_triangle_normals(std::vector<SbVec3f> &normals,
		     const struct rt_bot_internal *bot,
		     const std::vector<SbVec3f> &points,
		     const std::vector<Index> &triangles)
{
    if (!bot_authored_triangle_normals(normals, bot))
	generate_flat_triangle_normals(normals, points, triangles);
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
    bot_triangle_normals(normals, bot, points, indices);
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

static int
cad_mesh_part_geometry_from_bot(const struct rt_bot_internal *bot,
	Obol::PartGeometry &geometry)
{
    if (!bot || !bot->vertices || !bot->faces ||
	bot->num_vertices == 0 || bot->num_faces == 0)
	return 0;
    RT_BOT_CK_MAGIC(bot);

    Obol::TriMesh mesh;
    mesh.bounds.makeEmpty();
    mesh.positions.reserve(bot->num_vertices);
    for (size_t i = 0; i < bot->num_vertices; i++) {
	const SbVec3f point(static_cast<float>(bot->vertices[i * 3]),
	    static_cast<float>(bot->vertices[i * 3 + 1]),
	    static_cast<float>(bot->vertices[i * 3 + 2]));
	mesh.positions.push_back(point);
	mesh.bounds.extendBy(point);
    }

    if (bot->num_faces > std::numeric_limits<size_t>::max() / 3)
	return 0;
    mesh.indices.reserve(bot->num_faces * 3);
    for (size_t i = 0; i < bot->num_faces; i++) {
	const int *face = &bot->faces[i * 3];
	if (face[0] < 0 || face[1] < 0 || face[2] < 0 ||
	    static_cast<size_t>(face[0]) >= bot->num_vertices ||
	    static_cast<size_t>(face[1]) >= bot->num_vertices ||
	    static_cast<size_t>(face[2]) >= bot->num_vertices)
	    return 0;
	mesh.indices.push_back(static_cast<uint32_t>(face[0]));
	if (bot->orientation == RT_BOT_CW) {
	    mesh.indices.push_back(static_cast<uint32_t>(face[2]));
	    mesh.indices.push_back(static_cast<uint32_t>(face[1]));
	} else {
	    mesh.indices.push_back(static_cast<uint32_t>(face[1]));
	    mesh.indices.push_back(static_cast<uint32_t>(face[2]));
	}
    }
    if (mesh.bounds.isEmpty() || mesh.indices.empty())
	return 0;
    bot_triangle_normals(mesh.normals, bot, mesh.positions, mesh.indices);
    if (!canonicalize_corner_normal_mesh(mesh, mesh.normals))
	return 0;
    geometry.shaded = std::move(mesh);
    geometry.shadedCullBackfaces =
	bot->mode == RT_BOT_SOLID &&
	bot->orientation != RT_BOT_UNORIENTED &&
	bot->num_vertices <= static_cast<size_t>(INT_MAX) &&
	bot->num_faces <= static_cast<size_t>(INT_MAX) &&
	bg_trimesh_solid2(static_cast<int>(bot->num_vertices),
	    static_cast<int>(bot->num_faces),
	    const_cast<fastf_t *>(bot->vertices),
	    const_cast<int *>(bot->faces), NULL) == 0;
    return 1;
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

static void
sanitize_triangle_normals(std::vector<SbVec3f> &normals,
			  const std::vector<SbVec3f> &points,
			  const std::vector<int32_t> &triangles)
{
    if (normals.empty()) {
	generate_flat_triangle_normals(normals, points, triangles);
	return;
    }
    if (normals.size() != triangles.size())
	return;

    for (size_t i = 0; i + 2 < triangles.size(); i += 3) {
	const int32_t ia = triangles[i];
	const int32_t ib = triangles[i + 1];
	const int32_t ic = triangles[i + 2];
	if (ia < 0 || ib < 0 || ic < 0 ||
	    static_cast<size_t>(ia) >= points.size() ||
	    static_cast<size_t>(ib) >= points.size() ||
	    static_cast<size_t>(ic) >= points.size())
	    continue;

	SbVec3f faceNormal =
	    (points[static_cast<size_t>(ib)] - points[static_cast<size_t>(ia)]).cross(
		points[static_cast<size_t>(ic)] - points[static_cast<size_t>(ia)]);
	if (faceNormal.length() > 0.0f)
	    faceNormal.normalize();
	else
	    continue;

	float normalDot = 0.0f;
	for (size_t j = 0; j < 3; j++) {
	    SbVec3f &n = normals[i + j];
	    if (n.length() > 0.0f) {
		n.normalize();
		normalDot += n.dot(faceNormal);
	}
}

	if (normalDot < 0.0f) {
	    for (size_t j = 0; j < 3; j++)
		normals[i + j].negate();
	}

	for (size_t j = 0; j < 3; j++) {
	    if (normals[i + j].length() <= 0.0f)
		normals[i + j] = faceNormal;
	}
    }
}

static uint32_t
corner_normal_float_bits(float value)
{
    uint32_t bits = 0;
    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

struct CornerNormalVertexKey {
    uint32_t position;
    uint32_t normalX;
    uint32_t normalY;
    uint32_t normalZ;

    bool operator<(const CornerNormalVertexKey &other) const
    {
	if (position != other.position)
	    return position < other.position;
	if (normalX != other.normalX)
	    return normalX < other.normalX;
	if (normalY != other.normalY)
	    return normalY < other.normalY;
	return normalZ < other.normalZ;
    }
};

/* Obol's indexed mesh format binds one normal to each vertex.  BRL-CAD mesh
 * producers keep normals per triangle corner, so split only corners whose
 * normal differs instead of losing smooth shading or duplicating every vertex. */
static int
canonicalize_corner_normal_mesh(Obol::TriMesh &mesh,
	const std::vector<SbVec3f> &cornerNormals)
{
    if (cornerNormals.empty()) {
	mesh.normals.clear();
	return 1;
    }
    if (cornerNormals.size() != mesh.indices.size())
	return 0;

    std::map<CornerNormalVertexKey, uint32_t> vertexByCorner;
    std::vector<SbVec3f> positions;
    std::vector<SbVec3f> normals;
    std::vector<uint32_t> indices;
    positions.reserve(mesh.indices.size());
    normals.reserve(mesh.indices.size());
    indices.reserve(mesh.indices.size());
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
	const uint32_t sourceIndex = mesh.indices[i];
	if (sourceIndex >= mesh.positions.size())
	    return 0;
	SbVec3f normal = cornerNormals[i];
	if (normal.sqrLength() > 0.0f)
	    normal.normalize();
	else
	    normal.setValue(0.0f, 0.0f, 1.0f);
	const CornerNormalVertexKey key = {
	    sourceIndex,
	    corner_normal_float_bits(normal[0]),
	    corner_normal_float_bits(normal[1]),
	    corner_normal_float_bits(normal[2])
	};
	auto found = vertexByCorner.find(key);
	if (found == vertexByCorner.end()) {
	    const uint32_t index = static_cast<uint32_t>(positions.size());
	    positions.push_back(mesh.positions[sourceIndex]);
	    normals.push_back(normal);
	    found = vertexByCorner.emplace(key, index).first;
	}
	indices.push_back(found->second);
    }

    mesh.positions = std::move(positions);
    mesh.normals = std::move(normals);
    mesh.indices = std::move(indices);
    return 1;
}

static int
cad_mesh_part_geometry_from_indexed_face_set(
	const struct rt_primitive_indexed_face_set *faceSet,
	Obol::PartGeometry &geometry)
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
    static const char contract[] =
	"BObol-display-asset:shaded-triangles:brep-indexed-face-set-v2:";
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
	const struct bn_tol *tol, Obol::PartGeometry &geometry)
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
    shape->setIndexedTriangles(points.data(),
			       static_cast<int>(points.size()),
			       triangles.data(),
			       static_cast<int>(triangles.size()),
			       normals.empty() ? NULL : normals.data(),
			       static_cast<int>(normals.size()));
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

    struct bv_view_info viewInfo;
    source_view_info(&viewInfo, source);

    struct BObolMeshLodInfo info = BOBOL_MESH_LOD_INFO_INIT;
    if (bobol_mesh_lod_load_view(lod, &viewInfo, 0) >= 0 &&
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
    Obol::PartGeometry geometry;
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
	const struct bn_tol *tol, Obol::PartGeometry &geometry)
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
cad_mesh_append_hidden_line_edges(Obol::PartGeometry &geometry)
{
    if (!geometry.shaded)
	return 0;
    const Obol::TriMesh &mesh = *geometry.shaded;
    Obol::WireRep wire;
    wire.bounds = mesh.bounds;
    std::unordered_set<uint64_t> edges;
    uint32_t segmentId = 0;
    auto addEdge = [&](uint32_t a, uint32_t b) {
	if (a >= mesh.positions.size() || b >= mesh.positions.size())
	    return;
	if (a > b)
	    std::swap(a, b);
	const uint64_t key = (static_cast<uint64_t>(a) << 32) | b;
	if (!edges.insert(key).second)
	    return;
	wire.segmentPoints.push_back(mesh.positions[a]);
	wire.segmentPoints.push_back(mesh.positions[b]);
	wire.segmentIds.push_back(segmentId++);
    };
    for (size_t i = 0; i + 2 < mesh.indices.size(); i += 3) {
	addEdge(mesh.indices[i], mesh.indices[i + 1]);
	addEdge(mesh.indices[i + 1], mesh.indices[i + 2]);
	addEdge(mesh.indices[i + 2], mesh.indices[i]);
    }
    if (wire.segmentPoints.empty())
	return 0;
    geometry.wire = std::move(wire);
    return 1;
}

static int
cad_mesh_part_geometry_from_internal(struct rt_db_internal *intern,
	const SoBRLDatabaseSource *source, Obol::PartGeometry &geometry)
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
		 * service loads the view-selected PoP level. */
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

    struct rt_eval_wireframe_opts opts = RT_EVAL_WIREFRAME_OPTS_INIT;
    const int ret = rt_eval_wireframe(&vhead, &vlfree, dbip, path.c_str(),
				      &tol, &ttol, &opts);
    if (ret != BRLCAD_OK) {
	RT_FREE_VLIST(&vlfree, &vhead);
	bg_vlist_cleanup(&vlfree);
	return NULL;
    }

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
	convert_vlist(points, commands, &vhead);
    RT_FREE_VLIST(&vlfree, &vhead);
    bg_vlist_cleanup(&vlfree);
    if (points.empty() || points.size() != commands.size() ||
	points.size() > static_cast<size_t>(INT_MAX))
	return NULL;

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->setLineSet(points.data(), commands.data(),
		      static_cast<int>(points.size()));
    return shape;
}

static SoBRLMeshShape *
mesh_from_evaluated_points_path(SoBRLDatabaseSource *source)
{
    struct db_i *dbip = source ? source->getDatabase() : NULL;
    if (!source || !dbip)
	return NULL;

    const std::string path = database_source_evaluated_path_string(source);
    if (path.empty())
	return NULL;

    struct rt_primitive_indexed_face_set faceSet;
    memset(&faceSet, 0, sizeof(faceSet));
    const int ret = bobol_evaluated_points_evaluate_path_face_set(
			dbip, path.c_str(), &faceSet);
    if (ret != BRLCAD_OK)
	return NULL;

    SoBRLMeshShape *shape = mesh_from_indexed_face_set(&faceSet, source);
    bobol_evaluated_points_face_set_free(&faceSet);
    return shape;
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

    SoBRLMeshShape *shape = mesh_from_evaluated_points_path(source);
    if (!shape) {
	source->realizationDiagnostic =
	    "evaluated-points provider produced no drawable geometry";
	return -1;
    }

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
    SoBRLVListShape *sharedVListShape = NULL;
    SoBRLMeshShape *sharedMeshShape = NULL;
    SbBox3f cacheBounds;
    std::string cacheKey = realize_geometry_cache_key(dp);
    source_lod_cache_key_append(cacheKey, source, cacheBounds,
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP);
    std::map<std::string, SoBRLVListShape *>::iterator foundVList =
	cache->sharedMeshVListGeometry.find(cacheKey);
    if (foundVList != cache->sharedMeshVListGeometry.end()) {
	sharedVListShape = foundVList->second;
    } else {
	std::map<std::string, SoBRLMeshShape *>::iterator foundMesh =
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
	    SoSeparator *leaf = new SoSeparator;
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

    SoSeparator *leaf = new SoSeparator;
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
    SoBRLVListShape *sharedVListShape = NULL;
    SoBRLMeshShape *sharedMeshShape = NULL;
    SbBox3f cacheBounds;
    std::string cacheKey = realize_geometry_cache_key(dp);
    source_lod_cache_key_append(cacheKey, source, cacheBounds,
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP);
    std::map<std::string, SoBRLVListShape *>::iterator foundVList =
	cache->sharedMeshVListGeometry.find(cacheKey);
    if (foundVList != cache->sharedMeshVListGeometry.end()) {
	sharedVListShape = foundVList->second;
    } else {
	std::map<std::string, SoBRLMeshShape *>::iterator foundMesh =
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
	    SoSeparator *leaf = new SoSeparator;
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
	    Obol::PartGeometry generated;
	    int generatedGeometry = 0;
	    BObolSourceMeshRequest sourceMeshRequest;
	    bool lodBacked = false;
	    if (internalType == ID_PNTS) {
		generatedGeometry = cad_points_part_geometry_from_pnts(
		    static_cast<const struct rt_pnts_internal *>(
			validInternal.local.idb_ptr), generated);
	    } else if (wireGeometry) {
		generatedGeometry = cad_wire_part_geometry_from_lod_realization_internal(
		    &validInternal.local, source, localBounds, generated);
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
		    cad_source_mesh_request_from_bot(sourceMeshRequest, bot) &&
		    cad_wire_part_geometry_from_aabb(sourceMeshRequest.bounds,
			generated)) {
		    generatedGeometry = 1;
		    lodBacked = true;
		} else if (!generatedGeometry) {
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
			localBounds.isEmpty() ? NULL : &localBounds);
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
	if (!occurrence.geometry && sharedVListShape) {
	    Obol::PartGeometry generated;
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
	    Obol::PartGeometry generated;
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
	    occurrence.summary.lodActiveLevel = meshSummary.lodActiveLevel;
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
    return source->setCompactOccurrence(occurrence) > 0 ? 1 : 0;
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
	Obol::PartGeometry generated;
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
		Obol::PartGeometry generated;
		const int generatedGeometry = internalType == ID_PNTS ?
		    cad_points_part_geometry_from_pnts(
			static_cast<const struct rt_pnts_internal *>(
			    localIntern->idb_ptr), generated) :
		    (cad_wire_part_geometry_from_lod_realization_internal(
			localIntern, data->source, localBounds, generated) ||
		     cad_wire_part_geometry_from_plot_internal(localIntern,
			data->source, generated));
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
			localBounds.isEmpty() ? NULL : &localBounds);
		    cachedWire = data->cache->findMeshVListCadGeometry(cacheKey);
		}
	    } else if (data->source->lodBotThreshold.getValue() > 0) {
		Obol::PartGeometry generated;
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
			provisional.geometry =
			    std::make_shared<Obol::PartGeometry>(generated);
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
				data->source, tsp, provisionalPath.c_str(),
				dp->d_namep, primitive_type_label(localIntern),
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
	    input.occurrence.localTransform = mat_to_sbmatrix(tsp->ts_mat);
	    input.occurrence.summary = compact_occurrence_tree_summary(
		data->source, tsp, path, dp->d_namep,
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
	    if (!bounds.isEmpty()) {
		data->compact_bounds.extendBy(bounds);
		data->compact_bounds_valid = TRUE;
	    }
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
    std::map<std::string, SoBRLVListShape *>::iterator foundVList =
	data->cache->sharedMeshVListGeometry.find(cacheKey);
    if (foundVList != data->cache->sharedMeshVListGeometry.end()) {
	sharedVListShape = foundVList->second;
    } else {
	std::map<std::string, SoBRLMeshShape *>::iterator foundMesh =
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
	    if (!input.occurrence.geometry && sharedVListShape) {
		Obol::PartGeometry generated;
		if (cad_vlist_part_geometry(sharedVListShape, generated))
		    input.occurrence.geometry =
			data->cache->storeMeshVListCadGeometry(cacheKey,
			    std::move(generated), typeLabel, geometryKind);
	    }
	    input.occurrence.summary = compact_occurrence_tree_summary(
		data->source, tsp, path, dp->d_namep,
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
		Obol::PartGeometry generated;
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
		data->source, tsp, path, dp->d_namep, typeLabel, geometryKind,
		data->revision, BObolRealizedShapeSummary::SHAPE_MESH,
		static_cast<BObolMaterialColorSweep *>(data->material_sweep));
	    if (sharedMeshShape) {
		BObolRealizedShapeSummary meshSummary;
		realized_mesh_shape_summary(sharedMeshShape, meshSummary);
		input.occurrence.summary.lodAvailable = meshSummary.lodAvailable;
		input.occurrence.summary.lodActiveLevel =
		    meshSummary.lodActiveLevel;
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
	if (!bounds.isEmpty()) {
	    data->compact_bounds.extendBy(bounds);
	    data->compact_bounds_valid = TRUE;
	}
    } else {
	SoSeparator *leaf = realize_instance_leaf_separator(tsp);
	if (sharedVListShape) {
	    SoBRLVListShape *vlistShape = new SoBRLVListShape;
	    assign_realized_identity(vlistShape, tsp, path, dp->d_namep,
		typeLabel, data->revision, data->source);
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
		data->revision, data->source);
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
update_source_bounds_from_realized_geometry(SoBRLDatabaseSource *source)
{
    if (!source)
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
					   bounds.getMax());
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

static uint64_t
compact_next_revision(uint64_t revision)
{
    return revision == UINT64_MAX ? 1 : revision + 1;
}

static bool
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

static Obol::InstanceStyle
compact_effective_style(const BObolCompactInstanceEntry &entry)
{
    if (entry.highlighted)
	return entry.highlightedStyle;
    if (entry.selected)
	return entry.selectedStyle;
    return entry.normalStyle;
}

static void
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
    summary.visible = entry.visible;
    summary.selectable = entry.selectable;
    summary.selected = entry.selected;
    summary.highlighted = entry.highlighted;
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

static int
cad_part_key_for_geometry(const char *kind,
			  const Obol::PartGeometry &geometry,
			  std::string &key)
{
    if (!kind)
	return 0;

    CadGeometryHash hash;
    hash.appendString(kind);
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
			Obol::PartGeometry &geometry)
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
		       Obol::PartGeometry &geometry)
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
    return 1;
}

static std::string
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

static SbMatrix
cad_instance_matrix(const SoBRLDatabaseSource *source,
		    const SbMatrix &localMatrix)
{
    SbMatrix matrix = localMatrix;
    if (source && source->drawMatrixValid.getValue())
	matrix.multRight(source->drawMatrix.getValue());
    return matrix;
}

static SbMatrix
compact_mesh_asset_matrix(const SoBRLDatabaseSource *source,
	const BObolCompactInstanceEntry &entry)
{
    SbMatrix matrix = entry.sourceMeshRequest.meshAssetTransform;
    matrix.multRight(entry.placementTransform);
    return cad_instance_matrix(source, matrix);
}

struct cad_build_data {
    SoBRLDatabaseSource *source;
    SoBRLCadAssembly *assembly;
    std::map<std::string, Obol::PartId> partIdByKey;
    std::vector<Obol::PartUpdate> parts;
    std::vector<Obol::InstanceUpdate> instances;
    std::vector<Obol::InstanceId> hiddenInstances;
    std::vector<Obol::InstanceId> selectedInstances;
    std::vector<Obol::InstanceId> unpickableInstances;
    int ordinal;
    int unsupported;
    int wireCount;
    int shadedCount;
    int storeSemantics;
};

static void
cad_add_part_if_needed(cad_build_data &data,
		       const std::string &partKey,
		       const Obol::PartGeometry &geometry,
		       Obol::PartId &partId)
{
    std::map<std::string, Obol::PartId>::iterator found =
	data.partIdByKey.find(partKey);
    if (found != data.partIdByKey.end()) {
	partId = found->second;
	return;
    }

    partId = Obol::CadIdBuilder::hash128(partKey);
    data.partIdByKey[partKey] = partId;
    Obol::PartUpdate update;
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

static Obol::InstanceStyle
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

static const char *
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

static Obol::InstanceId
cad_source_parent_instance(const SoBRLDatabaseSource *source)
{
    if (!source)
	return Obol::CadIdBuilder::Root();
    const char *key = source->parentInstanceKey.getValue().getString();
    if (!key || !key[0])
	return Obol::CadIdBuilder::Root();
    return Obol::CadIdBuilder::hash128(key);
}

static uint8_t
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

static SoBRLCadAssembly::InstanceSemantic
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

static int
cad_proxy_corners(const BObolLodProxy &proxy, SbVec3f corners[8])
{
    if (!proxy.isValid())
	return 0;

    if (proxy.kind == BOBOL_LOD_PROXY_AABB) {
	const SbVec3f bmin = proxy.bounds.getMin();
	const SbVec3f bmax = proxy.bounds.getMax();
	corners[0].setValue(bmin[0], bmin[1], bmin[2]);
	corners[1].setValue(bmax[0], bmin[1], bmin[2]);
	corners[2].setValue(bmax[0], bmax[1], bmin[2]);
	corners[3].setValue(bmin[0], bmax[1], bmin[2]);
	corners[4].setValue(bmin[0], bmin[1], bmax[2]);
	corners[5].setValue(bmax[0], bmin[1], bmax[2]);
	corners[6].setValue(bmax[0], bmax[1], bmax[2]);
	corners[7].setValue(bmin[0], bmax[1], bmax[2]);
	return 1;
    }

    if (proxy.kind == BOBOL_LOD_PROXY_OBB) {
	SbVec3f ax = proxy.axisX;
	SbVec3f ay = proxy.axisY;
	SbVec3f az = proxy.axisZ;
	if (ax.length() > 0.0f)
	    ax.normalize();
	if (ay.length() > 0.0f)
	    ay.normalize();
	if (az.length() > 0.0f)
	    az.normalize();
	ax *= proxy.halfExtents[0];
	ay *= proxy.halfExtents[1];
	az *= proxy.halfExtents[2];
	corners[0] = proxy.center - ax - ay - az;
	corners[1] = proxy.center + ax - ay - az;
	corners[2] = proxy.center + ax + ay - az;
	corners[3] = proxy.center - ax + ay - az;
	corners[4] = proxy.center - ax - ay + az;
	corners[5] = proxy.center + ax - ay + az;
	corners[6] = proxy.center + ax + ay + az;
	corners[7] = proxy.center - ax + ay + az;
	return 1;
    }

    return 0;
}

static int
cad_wire_geometry_from_corners(const SbVec3f corners[8],
			       Obol::PartGeometry &geometry)
{
    static const int edges[12][2] = {
	{0, 1}, {1, 2}, {2, 3}, {3, 0},
	{4, 5}, {5, 6}, {6, 7}, {7, 4},
	{0, 4}, {1, 5}, {2, 6}, {3, 7}
    };

    Obol::WireRep wire;
    wire.bounds.makeEmpty();
    wire.segmentPoints.reserve(24);
    wire.segmentIds.reserve(12);
    for (size_t i = 0; i < 12; i++) {
	const SbVec3f &a = corners[edges[i][0]];
	const SbVec3f &b = corners[edges[i][1]];
	wire.segmentPoints.push_back(a);
	wire.segmentPoints.push_back(b);
	wire.segmentIds.push_back(static_cast<uint32_t>(i));
	wire.bounds.extendBy(a);
	wire.bounds.extendBy(b);
    }
    if (wire.bounds.isEmpty())
	return 0;
    geometry.wire = wire;
    return 1;
}

static int
cad_proxy_part_geometry(const BObolLodProxy &proxy,
			Obol::PartGeometry &geometry)
{
    SbVec3f corners[8];
    if (!cad_proxy_corners(proxy, corners))
	return 0;
    if (!cad_wire_geometry_from_corners(corners, geometry))
	return 0;
    /* This is a view-LoD AABB/OBB, not authored wire geometry.  Obol may
     * collapse it into a view-local point when the entire proxy is subpixel. */
    geometry.subpixelProxyEligible = true;
    geometry.structuralProxy = true;
    return 1;
}

static int
cad_mesh_payload_part_geometry(const BObolLodMeshPayload &payload,
			       int wire,
			       int shaded,
			       int shadedCullBackfaces,
			       Obol::PartGeometry &geometry,
			       BObolViewLodState::NormalStyle normalStyle =
				   BObolViewLodState::NORMAL_AUTHORED,
			       float normalCreaseAngle = 60.0f)
{
    if (!payload.isValid() || (!wire && !shaded))
	return 0;

    SbBox3f bounds;
    bounds.makeEmpty();
    for (size_t i = 0; i < payload.points.size(); i++)
	bounds.extendBy(payload.points[i]);
    if (bounds.isEmpty())
	return 0;

    if (shaded) {
	Obol::TriMesh mesh;
	mesh.positions = payload.points;
	mesh.bounds = bounds;
	mesh.indices.reserve(payload.coordIndex.size());
	for (size_t i = 0; i < payload.coordIndex.size(); i++) {
	    const int idx = payload.coordIndex[i];
	    if (idx < 0 || static_cast<size_t>(idx) >= payload.points.size())
		return 0;
	    mesh.indices.push_back(static_cast<uint32_t>(idx));
	}
	if (mesh.indices.empty())
	    return 0;
	std::vector<int32_t> cornerIndices;
	cornerIndices.reserve(payload.coordIndex.size());
	for (const int index : payload.coordIndex)
	    cornerIndices.push_back(index);
	std::vector<SbVec3f> cornerNormals;
	if (normalStyle == BObolViewLodState::NORMAL_FLAT) {
	    generate_flat_triangle_normals(cornerNormals, mesh.positions,
		cornerIndices);
	} else if (normalStyle == BObolViewLodState::NORMAL_SMOOTH) {
	    generate_smooth_triangle_normals(cornerNormals, mesh.positions,
		cornerIndices, normalCreaseAngle);
	} else if (payload.normals.size() == payload.coordIndex.size()) {
	    cornerNormals = payload.normals;
	}
	/* An indexed scan mesh without authored normals must stay indexed.
	 * Expanding flat corner normals here can turn N shared vertices into
	 * three vertices per face (tens of millions for Lucy), defeating PoP
	 * residency and creating multi-gigabyte transient maps.  Obol derives a
	 * flat face normal in the fragment shader, or per triangle in its
	 * immediate fallback, when this array remains empty. */
	if (!cornerNormals.empty()) {
	    sanitize_triangle_normals(cornerNormals, mesh.positions,
		cornerIndices);
	    if (!canonicalize_corner_normal_mesh(mesh, cornerNormals))
		return 0;
	}
	geometry.shaded = mesh;
	geometry.shadedCullBackfaces = shadedCullBackfaces != 0;
    }

    if (wire) {
	Obol::WireRep wireRep;
	wireRep.bounds = bounds;
	uint32_t edgeId = 0;
	wireRep.segmentPoints.reserve(payload.coordIndex.size() * 2);
	wireRep.segmentIds.reserve(payload.coordIndex.size());
	for (size_t i = 0; i + 2 < payload.coordIndex.size(); i += 3) {
	    int tri[3] = {
		payload.coordIndex[i],
		payload.coordIndex[i + 1],
		payload.coordIndex[i + 2]
	    };
	    for (int e = 0; e < 3; e++) {
		int a = tri[e];
		int b = tri[(e + 1) % 3];
		if (a < 0 || b < 0 ||
		    static_cast<size_t>(a) >= payload.points.size() ||
		    static_cast<size_t>(b) >= payload.points.size())
		    return 0;
		wireRep.segmentPoints.push_back(
		    payload.points[static_cast<size_t>(a)]);
		wireRep.segmentPoints.push_back(
		    payload.points[static_cast<size_t>(b)]);
		wireRep.segmentIds.push_back(edgeId++);
	    }
	}
	if (wireRep.segmentPoints.empty())
	    return 0;
	geometry.wire = wireRep;
    }

    /*
     * This is view-LoD presentation data with conservative full-mesh bounds.
     * When that entire extent is below one pixel, Obol may keep the retained
     * mesh and picking identity while drawing it through the aggregate point
     * channel.  Promotion is then a camera-local mask change, not a reload or
     * scene traversal.
     */
    geometry.subpixelProxyEligible = true;
    return 1;
}

static int
cad_progressive_mesh_part_geometry(
    const BObolLodProgressiveMeshPtr &progressive,
    int wire,
    int shaded,
    int shadedCullBackfaces,
    Obol::PartGeometry &geometry,
    BObolViewLodState::NormalStyle normalStyle =
	BObolViewLodState::NORMAL_AUTHORED,
    float normalCreaseAngle = 60.0f)
{
    if (!progressive || !progressive->isValid())
	return 0;
    const int residentLevel = progressive->residentLevel();
    /* residentLevel is the highest population prefix loaded from storage,
     * not necessarily the finest coordinate cut the retained exact arrays
     * can draw.  Once all points and faces needed by later PoP levels are
     * resident, those levels differ only in the renderer's coordinate snap
     * and require no further I/O.  Advertise that drawable frontier to Obol;
     * otherwise it clamps a correctly requested level (for example 9) back
     * to the last population-changing level (for example 4), leaving small
     * meshes visibly blocky forever. */
    int drawableLevel = residentLevel;
    for (int level = residentLevel + 1;
	 level <= progressive->maximumLevel(); ++level) {
	if (!progressive->canDrawLevel(level))
	    break;
	drawableLevel = level;
    }
    BObolLodMeshPayload payload;
    if (!progressive->copyLevel(payload, residentLevel) ||
	!cad_mesh_payload_part_geometry(payload, wire, shaded,
	    shadedCullBackfaces, geometry, normalStyle, normalCreaseAngle))
	return 0;

    const int minimumLevel = progressive->minimumLevel();
    const SbVec3f quantizationMinimum =
	progressive->quantizationMinimum();
    const SbVec3f quantizationMaximum =
	progressive->quantizationMaximum();
    const SbBox3f conservativeBounds(
	quantizationMinimum, quantizationMaximum);
    if (geometry.shaded.has_value()) {
	Obol::TriMesh &mesh = *geometry.shaded;
	/*
	 * A resident PoP prefix does not necessarily touch the asset's full
	 * extent.  Its quantization domain does, and is therefore the stable
	 * conservative bound required by view culling and the aggregate
	 * subpixel proxy.  Prefix-derived bounds caused objects to cross the
	 * proxy threshold as data arrived, repeatedly rebuilding large retained
	 * submissions at an unchanged camera.
	 */
	if (!conservativeBounds.isEmpty())
	    mesh.bounds = conservativeBounds;
	mesh.progressiveMinimumLevel =
	    static_cast<uint8_t>(std::max(0, minimumLevel));
	mesh.progressiveResidentLevel =
	    static_cast<uint8_t>(std::max(0, drawableLevel));
	mesh.progressiveQuantizationMinimum = quantizationMinimum;
	mesh.progressiveQuantizationMaximum = quantizationMaximum;
	for (int level = 0; level < BOBOL_MESH_LOD_LEVEL_COUNT; ++level) {
	    const size_t count = progressive->faceCount(level) * 3;
	    mesh.progressiveIndexCount[level] = static_cast<uint32_t>(
		std::min(count,
		    static_cast<size_t>(std::numeric_limits<uint32_t>::max())));
	}
	/* Record the vertex extent of each cumulative index prefix once while
	 * constructing the geometry.  Recomputing it from every visible part on
	 * every render is quadratic bookkeeping for large distinct-mesh scenes.
	 * Derive it from the finalized Obol indices rather than the source point
	 * counts because normal canonicalization may split shared vertices. */
	size_t scannedIndexCount = 0;
	uint32_t maximumIndex = 0;
	for (int level = 0; level < BOBOL_MESH_LOD_LEVEL_COUNT; ++level) {
	    const size_t indexCount = std::min<size_t>(
		mesh.progressiveIndexCount[level], mesh.indices.size());
	    for (; scannedIndexCount < indexCount; ++scannedIndexCount)
		maximumIndex = std::max(maximumIndex,
		    mesh.indices[scannedIndexCount]);
	    const size_t positionCount = indexCount ?
		static_cast<size_t>(maximumIndex) + 1 : 0;
	    mesh.progressivePositionCount[level] = static_cast<uint32_t>(
		std::min(positionCount,
		    static_cast<size_t>(std::numeric_limits<uint32_t>::max())));
	}
    }
    if (geometry.wire.has_value()) {
	Obol::WireRep &wireRep = *geometry.wire;
	if (!conservativeBounds.isEmpty())
	    wireRep.bounds = conservativeBounds;
	wireRep.progressiveMinimumLevel =
	    static_cast<uint8_t>(std::max(0, minimumLevel));
	wireRep.progressiveResidentLevel =
	    static_cast<uint8_t>(std::max(0, drawableLevel));
	wireRep.progressiveQuantizationMinimum = quantizationMinimum;
	wireRep.progressiveQuantizationMaximum = quantizationMaximum;
	for (int level = 0; level < BOBOL_MESH_LOD_LEVEL_COUNT; ++level) {
	    const size_t count = progressive->faceCount(level) * 3;
	    wireRep.progressiveSegmentCount[level] = static_cast<uint32_t>(
		std::min(count,
		    static_cast<size_t>(std::numeric_limits<uint32_t>::max())));
	}
    }
    return 1;
}

static bool
cad_progressive_geometry_can_draw(const Obol::PartGeometry *geometry,
				  int wire, int shaded, int activeLevel)
{
    if (!geometry || activeLevel < 0)
	return false;
    if (wire) {
	if (!geometry->wire.has_value() ||
	    !geometry->wire->isProgressive() ||
	    geometry->wire->progressiveResidentLevel <
		static_cast<uint8_t>(std::min(15, activeLevel)))
	    return false;
    }
    if (shaded) {
	if (!geometry->shaded.has_value() ||
	    !geometry->shaded->isProgressive() ||
	    geometry->shaded->progressiveResidentLevel <
		static_cast<uint8_t>(std::min(15, activeLevel)))
	    return false;
    }
    return wire || shaded;
}

static std::string
cad_view_lod_assembly_key(const SoBRLDatabaseSource *source,
			  const BObolViewLodState::CadPayload *payload,
			  const SbMatrix &matrix)
{
    std::string key = source ? source->path.getValue().getString() : "";
    key += "|";
    key += source ? source->instanceKey.getValue().getString() : "";
    key += "|";
    key += payload ? payload->cacheKey.getString() : "";
    char buf[256] = {0};
    snprintf(buf, sizeof(buf), "|%d|%d|%llu|%llu|%u|%u|%d|%d|%d|%d|%d|%d|%d|%g|%d",
	     payload ? payload->resultKind : 0,
	     payload ? payload->qualityTier : 0,
	     payload ? (unsigned long long)payload->viewRevision : 0ULL,
	     payload ? (unsigned long long)payload->policyRevision : 0ULL,
	     source ? source->sourceRevision.getValue() : 0,
	     source ? source->inputsRevision.getValue() : 0,
	     source ? source->drawMode.getValue() : 0,
	     source ? source->visible.getValue() : 0,
	     source ? source->selected.getValue() : 0,
	     source ? source->highlighted.getValue() : 0,
	     source ? source->colorOverride.getValue() : 0,
	     source ? source->materialColorValid.getValue() : 0,
	     source ? source->lineWidth.getValue() : 0,
	     source ? source->transparency.getValue() : 0.0f,
	     source ? source->drawMatrixValid.getValue() : 0);
    key += buf;
    const SbColor selectedColor = source ?
	source->selectedColor.getValue() : SbColor(1.0f, 1.0f, 1.0f);
    const SbColor highlightedColor = source ?
	source->highlightedColor.getValue() : SbColor(1.0f, 1.0f, 0.0f);
    const SbColor ghostedColor = source ?
	source->ghostedColor.getValue() : SbColor(0.55f, 0.55f, 0.55f);
    for (int i = 0; i < 3; i++) {
	char cbuf[128] = {0};
	snprintf(cbuf, sizeof(cbuf), "|sc=%.9g|hc=%.9g|gc=%.9g",
		 selectedColor[i], highlightedColor[i], ghostedColor[i]);
	key += cbuf;
    }
    const float *m = matrix[0];
    for (int i = 0; i < 16; i++) {
	char mbuf[64] = {0};
	snprintf(mbuf, sizeof(mbuf), "|%.9g", m[i]);
	key += mbuf;
    }
    return key;
}

static SoBRLCadAssembly *
cad_view_lod_assembly(const SoBRLDatabaseSource *source,
		      const BObolViewLodState::CadPayload *payload,
		      const BObolViewLodState *viewState)
{
    if (!source || !payload || !viewState || !payload->isValid() ||
	!source->visible.getValue())
	return NULL;

    const SbMatrix matrix = cad_instance_matrix(source, SbMatrix::identity());
    const BObolViewLodState::NormalStyle normalStyle =
	viewState->getNormalStyle();
    const float normalCreaseAngle = viewState->getNormalCreaseAngle();
    std::string assemblyKey =
	cad_view_lod_assembly_key(source, payload, matrix);
    assemblyKey += "|normal-style=";
    assemblyKey += std::to_string(static_cast<int>(normalStyle));
    assemblyKey += "|normal-crease=";
    assemblyKey += std::to_string(normalCreaseAngle);
    SbString cachedKey;
    SoBRLCadAssembly *assembly = dynamic_cast<SoBRLCadAssembly *>(
	viewState->findCadPresentation(source, &cachedKey));
    if (assembly &&
	bu_strcmp(cachedKey.getString(), assemblyKey.c_str()) == 0)
	return assembly;

    Obol::PartGeometry geometry;
    std::shared_ptr<const Obol::PartGeometry> preparedGeometry;
    int shadedCount = 0;
    int wireCount = 0;
    SoBRLPickDetail::PrimitiveKind primitiveKind =
	SoBRLPickDetail::LINE_SEGMENT;
    const int sourceDrawMode = source_record_draw_mode(source);
    if (payload->resultKind == BOBOL_LOD_RESULT_AABB ||
	payload->resultKind == BOBOL_LOD_RESULT_PROXY) {
	if (!cad_proxy_part_geometry(payload->proxy, geometry))
	    return NULL;
	wireCount = 1;
    } else {
	const bool hiddenLine =
	    sourceDrawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
	const int wire =
	    (sourceDrawMode == BOBOL_LOD_DRAW_WIRE ||
	     hiddenLine ||
	     payload->drawMode == BOBOL_LOD_DRAW_WIRE) ? 1 : 0;
	const int shaded = hiddenLine || !wire ? 1 : 0;
	if (payload->progressiveMesh &&
	    normalStyle == BObolViewLodState::NORMAL_AUTHORED &&
	    payload->preparedCadGeometry &&
	    payload->preparedCadGeometryRevision ==
		payload->progressiveMesh->revision() &&
	    cad_progressive_geometry_can_draw(
		payload->preparedCadGeometry.get(), wire, shaded,
		payload->activeLevel))
	    preparedGeometry = payload->preparedCadGeometry;
	if (!preparedGeometry &&
	    !cad_mesh_payload_part_geometry(payload->mesh, wire, shaded,
		payload->shadedCullBackfaces,
		geometry, normalStyle, normalCreaseAngle))
	    return NULL;
	const Obol::PartGeometry *publishedGeometry =
	    preparedGeometry ? preparedGeometry.get() : &geometry;
	if (publishedGeometry->wire.has_value())
	    wireCount = 1;
	if (publishedGeometry->shaded.has_value()) {
	    shadedCount = 1;
	    primitiveKind = SoBRLPickDetail::FACE;
	}
    }

    if (!assembly) {
	assembly = new SoBRLCadAssembly;
	viewState->setCadPresentation(source, assembly);
    }
    assembly->beginUpdate();
    assembly->clear();
    assembly->clearSemanticMap();

    std::string partKey = "view-lod:";
    partKey += source->path.getValue().getString();
    partKey += ":";
    partKey += payload->cacheKey.getString();
    partKey += ":";
    partKey += payload->resultKind == BOBOL_LOD_RESULT_AABB ? "aabb" :
	       (payload->resultKind == BOBOL_LOD_RESULT_PROXY ? "proxy" :
		(wireCount ? "wire-mesh" : "mesh"));
    Obol::PartId partId = Obol::CadIdBuilder::hash128(partKey);
    Obol::SharedPartUpdate part;
    part.part = partId;
    part.geometry = preparedGeometry ? preparedGeometry :
	std::make_shared<const Obol::PartGeometry>(std::move(geometry));
    assembly->upsertSharedParts(
	std::vector<Obol::SharedPartUpdate>(1, std::move(part)));

    std::string instanceKey =
	cad_instance_key(source, source->path.getValue().getString(), 0);
    instanceKey += ":view-lod:";
    instanceKey += payload->cacheKey.getString();
    Obol::InstanceId instanceId = Obol::CadIdBuilder::hash128(instanceKey);
    Obol::InstanceRecord record;
    record.part = partId;
    record.localToRoot = matrix;
    record.parent = cad_source_parent_instance(source);
    record.childName = cad_source_leaf_name(source);
    record.occurrenceIndex = source->occurrenceIndex.getValue();
    record.boolOp = cad_source_boolean_operation(source);
    record.style = cad_source_style(source);
    Obol::InstanceUpdate update;
    update.instance = instanceId;
    update.record = record;
    assembly->upsertInstances(std::vector<Obol::InstanceUpdate>(1, update));
    assembly->setInstanceSemantic(instanceId,
	cad_source_semantic(source, primitiveKind));

    if (source_record_draw_mode(source) == BOBOL_LOD_DRAW_HIDDEN_LINE)
	assembly->drawMode = SoCADAssembly::HIDDEN_LINE;
    else if (shadedCount > 0 && wireCount > 0)
	assembly->drawMode = SoCADAssembly::SHADED_WITH_EDGES;
    else if (shadedCount > 0)
	assembly->drawMode = SoCADAssembly::SHADED;
    else
	assembly->drawMode = SoCADAssembly::WIREFRAME;
    if (payload->resultKind == BOBOL_LOD_RESULT_AABB ||
	payload->resultKind == BOBOL_LOD_RESULT_PROXY)
	assembly->pickMode = SoCADAssembly::PICK_BOUNDS;
    else
	assembly->pickMode = SoCADAssembly::PICK_AUTO;
    assembly->endUpdate();

    viewState->setCadPresentation(source, assembly, assemblyKey.c_str());
    return assembly;
}

static bool
cad_lod_path_equal(const SbString &a, const SbString &b)
{
    const char *ap = a.getString();
    const char *bp = b.getString();
    if (!ap || !bp)
	return false;
    if (ap[0] == '/')
	ap++;
    if (bp[0] == '/')
	bp++;
    return bu_strcmp(ap, bp) == 0;
}

static const BObolViewLodState::CadPayload *
cad_compact_payload_for_entry(
    const std::vector<const BObolViewLodState::CadPayload *> &payloads,
    const BObolCompactInstanceEntry &entry)
{
    const SbString &occurrenceKey = compact_instance_identity(entry);
    for (const BObolViewLodState::CadPayload *payload : payloads) {
	if (payload && payload->sourceInstanceKey.getLength() > 0 &&
	    bu_strcmp(payload->sourceInstanceKey.getString(),
		occurrenceKey.getString()) == 0)
	    return payload;
    }
    /* Source-wide legacy results deliberately have no occurrence key.  They
     * may still bind by path, but an occurrence-specific result never aliases
     * a sibling merely because both paths are textual matches. */
    for (const BObolViewLodState::CadPayload *payload : payloads) {
	if (payload && payload->sourceInstanceKey.getLength() == 0 &&
	    cad_lod_path_equal(payload->sourcePath, entry.semantic.path))
	    return payload;
    }
    return NULL;
}

SoBRLCadAssembly *
SoBRLDatabaseSource::currentCompactViewLodAssembly(
    const BObolViewLodState *viewState) const
{
    if (!this->d->compactIndex || !viewState)
	return NULL;
    SoBRLCadAssembly *assembly = dynamic_cast<SoBRLCadAssembly *>(
	viewState->findCadPresentation(this));
    if (!assembly || !assembly->compactPresentationInitialized ||
	assembly->compactPresentationIndex != this->d->compactIndex ||
	assembly->compactPresentationSourceRevision !=
	    this->sourceRevision.getValue() ||
	assembly->compactPresentationInputsRevision !=
	    this->inputsRevision.getValue() ||
	assembly->compactPresentationCadBatchRevision !=
	    this->cadBatchRevisionGet() ||
	assembly->compactPresentationDrawMode != source_record_draw_mode(this) ||
	assembly->compactPresentationPayloadRevision !=
	    viewState->cadRevision())
	return NULL;
    return assembly;
}

int
SoBRLDatabaseSource::getCompactViewLodSupersededFallbackCount(
    const BObolViewLodState *viewState, std::vector<SbString> *paths) const
{
    if (paths)
	paths->clear();
    if (!this->d->compactIndex || !viewState)
	return 0;

    SoBRLCadAssembly *assembly = dynamic_cast<SoBRLCadAssembly *>(
	viewState->findCadPresentation(this));
    if (!assembly || !assembly->compactPresentationInitialized ||
	assembly->compactPresentationIndex != this->d->compactIndex)
	return 0;

    std::vector<const BObolViewLodState::CadPayload *> payloads;
    viewState->findCadPayloads(this, payloads);
    std::unordered_set<std::string_view> payloadKeys;
    payloadKeys.reserve(payloads.size());
    for (const BObolViewLodState::CadPayload *payload : payloads)
	if (payload && payload->isValid() &&
	    payload->sourceInstanceKey.getLength() > 0)
	    payloadKeys.emplace(payload->sourceInstanceKey.getString(),
		static_cast<size_t>(
		    payload->sourceInstanceKey.getLength()));

    int count = 0;
    for (const BObolCompactInstanceEntry &entry :
	 this->d->compactIndex->entries) {
	const char *kind = entry.shapeSummary.geometryKind.getString();
	if (!BU_STR_EQUAL(kind, "aabb") && !BU_STR_EQUAL(kind, "obb"))
	    continue;
	const SbString &occurrenceKey = compact_instance_identity(entry);
	const std::string_view occurrenceView(
	    occurrenceKey.getString(),
	    static_cast<size_t>(occurrenceKey.getLength()));
	if (payloadKeys.find(occurrenceView) == payloadKeys.end())
	    continue;
	const auto presentation =
	    assembly->compactInstancePresentations.find(entry.instance);
	if (presentation == assembly->compactInstancePresentations.end() ||
	    presentation->second.activePart != entry.part)
	    continue;
	count++;
	if (paths)
	    paths->push_back(entry.semantic.path);
    }
    return count;
}

int
SoBRLDatabaseSource::getCompactViewLodActiveFallbackCount(
    const BObolViewLodState *viewState, std::vector<SbString> *paths) const
{
    if (paths)
	paths->clear();
    if (!this->d->compactIndex || !viewState)
	return 0;

    SoBRLCadAssembly *assembly = dynamic_cast<SoBRLCadAssembly *>(
	viewState->findCadPresentation(this));
    if (!assembly || !assembly->compactPresentationInitialized ||
	assembly->compactPresentationIndex != this->d->compactIndex)
	return 0;

    int count = 0;
    for (const BObolCompactInstanceEntry &entry :
	 this->d->compactIndex->entries) {
	const char *kind = entry.shapeSummary.geometryKind.getString();
	if (!BU_STR_EQUAL(kind, "aabb") && !BU_STR_EQUAL(kind, "obb"))
	    continue;
	const auto presentation =
	    assembly->compactInstancePresentations.find(entry.instance);
	if (presentation == assembly->compactInstancePresentations.end() ||
	    presentation->second.activePart != entry.part)
	    continue;
	count++;
	if (paths)
	    paths->push_back(entry.semantic.path);
    }
    return count;
}

SoBRLCadAssembly *
SoBRLDatabaseSource::compactViewLodAssembly(
    const std::vector<const BObolViewLodState::CadPayload *> &payloads,
    const BObolViewLodState *viewState)
    const
{
    const bool presentationTiming =
	getenv("BOBOL_COMPACT_PRESENTATION_TIMING") != NULL;
    const int64_t presentationStarted =
	presentationTiming ? bu_gettime() : 0;
    if (!this->d->compactIndex || !viewState)
	return NULL;

    SoBRLCadAssembly *assembly = dynamic_cast<SoBRLCadAssembly *>(
	viewState->findCadPresentation(this));
    if (!assembly) {
	assembly = new SoBRLCadAssembly;
	viewState->setCadPresentation(this, assembly);
    }
    assembly->reserveCompactPresentationCapacity(std::max(
	this->d->compactExpectedInstanceCount,
	this->d->compactIndex->entries.size()));

    const int sourceDrawMode = source_record_draw_mode(this);
    const uint64_t payloadRevision = viewState->cadRevision();
    const BObolViewLodState::NormalStyle normalStyle =
	viewState->getNormalStyle();
    const float normalCreaseAngle = viewState->getNormalCreaseAngle();
    if (assembly->compactPresentationInitialized &&
	assembly->compactPresentationIndex == this->d->compactIndex &&
	assembly->compactPresentationSourceRevision ==
	    this->sourceRevision.getValue() &&
	assembly->compactPresentationInputsRevision ==
	    this->inputsRevision.getValue() &&
	assembly->compactPresentationCadBatchRevision ==
	    this->cadBatchRevisionGet() &&
	assembly->compactPresentationDrawMode == sourceDrawMode &&
	assembly->compactPresentationPayloadRevision == payloadRevision)
	return assembly;
    const bool hiddenLine = sourceDrawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    bool reset = !assembly->compactPresentationInitialized ||
	assembly->compactPresentationIndex != this->d->compactIndex ||
	assembly->compactPresentationSourceRevision !=
	    this->sourceRevision.getValue() ||
	assembly->compactPresentationInputsRevision !=
	    this->inputsRevision.getValue() ||
	assembly->compactPresentationDrawMode != sourceDrawMode;

    /*
     * Streaming a large root appends/replaces a bounded occurrence batch in
     * the same compact index.  Clearing and reinserting the whole compiled
     * assembly on every batch made population O(n^2), invalidated every part
     * generation, and forced the renderer to re-snap every prior PoP cut.
     * Preserve the assembly when every old occurrence still exists and route
     * only new or geometry-changed records through the normal upsert path.
     * A removal or index replacement remains an authoritative full reset.
     */
    std::vector<size_t> structuralEntryIndices;
    const bool cadBatchChanged = !reset &&
	assembly->compactPresentationCadBatchRevision !=
	    this->cadBatchRevisionGet();
    SbBool sparseCadBatchUpdate = FALSE;
    if (cadBatchChanged)
	sparseCadBatchUpdate = this->getCadBatchChangedEntries(
	    assembly->compactPresentationCadBatchRevision,
	    structuralEntryIndices);
    if (cadBatchChanged && !sparseCadBatchUpdate) {
	structuralEntryIndices.reserve(
	    this->d->compactIndex->entries.size() >
		assembly->compactInstancePresentations.size() ?
	    this->d->compactIndex->entries.size() -
		assembly->compactInstancePresentations.size() : 0);
	size_t retainedPresentations = 0;
	for (size_t i = 0; i < this->d->compactIndex->entries.size(); i++) {
	    const BObolCompactInstanceEntry &entry =
		this->d->compactIndex->entries[i];
	    const auto found =
		assembly->compactInstancePresentations.find(entry.instance);
	    if (found == assembly->compactInstancePresentations.end()) {
		structuralEntryIndices.push_back(i);
	    } else {
		retainedPresentations++;
		const SoBRLCadAssembly::CompactInstancePresentation &prior =
		    found->second;
		if (prior.geometryRevision != entry.geometryRevision ||
		    prior.appearanceRevision != entry.appearanceRevision ||
		    prior.placementRevision != entry.placementRevision ||
		    prior.visibilityRevision != entry.visibilityRevision ||
		    prior.selectionRevision != entry.selectionRevision)
		    structuralEntryIndices.push_back(i);
	    }
	}
	if (retainedPresentations !=
		assembly->compactInstancePresentations.size()) {
	    reset = true;
	    structuralEntryIndices.clear();
	}
    } else if (sparseCadBatchUpdate) {
	/*
	 * The source mutation supplied an exact, non-consuming occurrence
	 * journal.  Validate only those records; other occurrences are known to
	 * match the presentation acknowledged at the prior CAD-batch revision.
	 */
	for (auto it = structuralEntryIndices.begin();
	     it != structuralEntryIndices.end();) {
	    if (*it >= this->d->compactIndex->entries.size()) {
		reset = true;
		structuralEntryIndices.clear();
		break;
	    }
	    const BObolCompactInstanceEntry &entry =
		this->d->compactIndex->entries[*it];
	    const auto prior =
		assembly->compactInstancePresentations.find(entry.instance);
	    if (prior == assembly->compactInstancePresentations.end()) {
		++it;
		continue;
	    }
	    const SoBRLCadAssembly::CompactInstancePresentation &presented =
		prior->second;
	    if (presented.geometryRevision == entry.geometryRevision &&
		presented.appearanceRevision == entry.appearanceRevision &&
		presented.placementRevision == entry.placementRevision &&
		presented.visibilityRevision == entry.visibilityRevision &&
		presented.selectionRevision == entry.selectionRevision)
		it = structuralEntryIndices.erase(it);
	    else
		++it;
	}
    }

    /* A camera epoch usually changes only the active PoP level on a bounded
     * wave of occurrences.  Preserve the compiled assembly and route those
     * changes directly by occurrence key.  Source structure/style revisions
     * still take the authoritative full-scan path. */
    const bool payloadOnlyUpdate = !reset &&
	assembly->compactPresentationCadBatchRevision ==
	    this->cadBatchRevisionGet();
    const bool incrementalStructureUpdate =
	!reset && !payloadOnlyUpdate;
    std::vector<SbString> changedOccurrenceKeys;
    SbBool fullPayloadResync = TRUE;
    if (!reset)
	viewState->cadOccurrenceChangesSince(
	    this, assembly->compactPresentationPayloadRevision,
	    changedOccurrenceKeys, fullPayloadResync);
    const bool sparsePayloadUpdate = !reset && !fullPayloadResync;
    if (sparsePayloadUpdate && changedOccurrenceKeys.empty()) {
	if (incrementalStructureUpdate && !structuralEntryIndices.empty()) {
	    /* Continue below with just the compact-index delta. */
	} else {
	    if (incrementalStructureUpdate)
		assembly->compactPresentationCadBatchRevision =
		    this->cadBatchRevisionGet();
	    assembly->compactPresentationPayloadRevision = payloadRevision;
	    viewState->acknowledgeCadOccurrenceChanges(this, payloadRevision);
	    return assembly;
	}
    }

    std::vector<const BObolViewLodState::CadPayload *>
	authoritativePayloads;
    const std::vector<const BObolViewLodState::CadPayload *> *
	payloadInput = &payloads;
    if (!sparsePayloadUpdate && payloads.empty()) {
	viewState->findCadPayloadsUnordered(this, authoritativePayloads);
	payloadInput = &authoritativePayloads;
    }

    /* Both collections may contain thousands of records.  Route modern
     * occurrence-specific payloads by key instead of scanning every payload
     * for every compact leaf.  Retain the small path scan only for legacy
     * source-wide payloads with no occurrence identity. */
    std::unordered_map<std::string_view,
	const BObolViewLodState::CadPayload *> payloadByInstance;
    std::vector<const BObolViewLodState::CadPayload *> sourceWidePayloads;
    payloadByInstance.reserve(payloadInput->size());
    for (const BObolViewLodState::CadPayload *payload : *payloadInput) {
	/*
	 * Payload input comes from the validated view-state occurrence table or
	 * from the caller's result lookup.  It cannot become invalid in place:
	 * replacement and eviction publish/erase whole shared payload objects.
	 */
	if (!payload)
	    continue;
	if (payload->sourceInstanceKey.getLength() > 0)
	    payloadByInstance[std::string_view(
		payload->sourceInstanceKey.getString(),
		static_cast<size_t>(
		    payload->sourceInstanceKey.getLength()))] = payload;
	else
	    sourceWidePayloads.push_back(payload);
    }

    std::vector<size_t> entryIndices;
    if (sparsePayloadUpdate || incrementalStructureUpdate) {
	entryIndices = structuralEntryIndices;
	entryIndices.reserve(entryIndices.size() +
	    changedOccurrenceKeys.size());
	for (const SbString &occurrenceKey : changedOccurrenceKeys) {
	    const auto found =
		this->d->compactIndex->entryIndexByKey.find(
		    occurrenceKey.getString());
	    if (found != this->d->compactIndex->entryIndexByKey.end())
		entryIndices.push_back(found->second);
	}
	std::sort(entryIndices.begin(), entryIndices.end());
	entryIndices.erase(
	    std::unique(entryIndices.begin(), entryIndices.end()),
	entryIndices.end());
    } else {
	entryIndices.resize(this->d->compactIndex->entries.size());
	std::iota(entryIndices.begin(), entryIndices.end(), 0);
    }
    const int64_t presentationIndexed =
	presentationTiming ? bu_gettime() : 0;

    std::vector<Obol::SharedPartUpdate> baseParts;
    std::vector<Obol::SharedPartUpdate> lodSharedParts;
    std::vector<Obol::InstanceUpdate> instances;
    std::vector<Obol::InstanceStyleUpdate> instanceStyles;
    std::vector<Obol::InstanceLodUpdate> lodLevelUpdates;
    std::unordered_set<Obol::PartId, std::hash<Obol::PartId>>
	progressivePartsUpdated;
    std::unordered_set<Obol::PartId, std::hash<Obol::PartId>> partsToRemove;
    bool instanceSetsChanged = reset;

    if (reset) {
	baseParts.reserve(this->d->compactIndex->parts.size());
	assembly->compactPartChannels.clear();
	for (const BObolCompactPartReference &partRef :
	     this->d->compactIndex->parts) {
	    if (!partRef.geometry)
		continue;
	    Obol::SharedPartUpdate part;
	    part.part = partRef.part;
	    part.geometry = partRef.geometry;
	    assembly->compactPartChannels[part.part] =
		(part.geometry->wire.has_value() ? 1u : 0u) |
		(part.geometry->shaded.has_value() ? 2u : 0u) |
		(part.geometry->points.has_value() ? 4u : 0u);
	    baseParts.push_back(std::move(part));
	}
	instances = this->d->compactIndex->instances;
	assembly->compactInstancePresentations.clear();
	assembly->compactActivePartReferences.clear();
	assembly->compactLodParts.clear();
    } else if (!structuralEntryIndices.empty()) {
	std::unordered_set<Obol::PartId, std::hash<Obol::PartId>>
	    addedBaseParts;
	addedBaseParts.reserve(structuralEntryIndices.size());
	for (const size_t i : structuralEntryIndices) {
	    if (i >= this->d->compactIndex->entries.size())
		continue;
	    const BObolCompactInstanceEntry &entry =
		this->d->compactIndex->entries[i];
	    const auto prior =
		assembly->compactInstancePresentations.find(entry.instance);
	    if (prior != assembly->compactInstancePresentations.end() &&
		prior->second.geometryRevision == entry.geometryRevision)
		continue;
	    if (!entry.geometry || !addedBaseParts.insert(entry.part).second)
		continue;
	    Obol::SharedPartUpdate part;
	    part.part = entry.part;
	    part.geometry = entry.geometry;
	    baseParts.push_back(std::move(part));
	    assembly->compactPartChannels[entry.part] =
		(entry.geometry->wire.has_value() ? 1u : 0u) |
		(entry.geometry->shaded.has_value() ? 2u : 0u) |
		(entry.geometry->points.has_value() ? 4u : 0u);
	}
    }

    const bool incrementalUpdate = !reset;
    size_t wireCount = incrementalUpdate ?
	assembly->compactWirePresentationCount : 0;
    size_t shadedCount = incrementalUpdate ?
	assembly->compactShadedPresentationCount : 0;
    for (const size_t i : entryIndices) {
	const BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[i];
	const SbString &occurrenceKey = compact_instance_identity(entry);
	const BObolViewLodState::CadPayload *payload = NULL;
	if (sparsePayloadUpdate) {
	    payload = viewState->findCadForOccurrence(
		this, occurrenceKey);
	} else {
	    const auto payloadIt = payloadByInstance.find(
		std::string_view(occurrenceKey.getString(),
		    static_cast<size_t>(occurrenceKey.getLength())));
	    payload =
		payloadIt != payloadByInstance.end() ? payloadIt->second :
		cad_compact_payload_for_entry(sourceWidePayloads, entry);
	}

	const auto existingPresentation =
	    assembly->compactInstancePresentations.find(entry.instance);
	const bool presentationWasPresent =
	    existingPresentation !=
		assembly->compactInstancePresentations.end();
	SoBRLCadAssembly::CompactInstancePresentation &presentation =
	    presentationWasPresent ? existingPresentation->second :
	    assembly->compactInstancePresentations[entry.instance];
	const uint8_t previousChannels = presentation.channels;
	const Obol::PartId previousPart = presentation.activePart;
	if (reset) {
	    presentation.activePart = entry.part;
	    presentation.channels = assembly->compactPartChannels[entry.part];
	    presentation.activeLevel = -1;
	}

    std::string payloadKey;
	Obol::PartId desiredPart = entry.part;
	uint8_t desiredChannels = assembly->compactPartChannels[entry.part];
	bool desiredGeometryValid = false;
	int desiredActiveLevel = -1;
	if (payload) {
	    payloadKey.reserve(
		static_cast<size_t>(payload->cacheKey.getLength()) + 128u);
	    payloadKey.assign(payload->cacheKey.getString(),
		static_cast<size_t>(payload->cacheKey.getLength()));
	    payloadKey += ':';
	    payloadKey += std::to_string(payload->resultKind);
	    payloadKey += ':';
	    payloadKey += std::to_string(payload->qualityTier);
	    payloadKey += ':';
	    if (payload->progressiveMesh) {
		payloadKey += "resident-revision=";
		payloadKey += std::to_string(
		    payload->progressiveMesh->revision());
		desiredActiveLevel = payload->activeLevel;
	    } else {
		payloadKey += "level=";
		payloadKey += std::to_string(payload->activeLevel);
	    }
	    payloadKey += ':';
	    payloadKey += std::to_string(sourceDrawMode);
	    payloadKey += ":normal-style=";
	    payloadKey += std::to_string(static_cast<int>(normalStyle));
	    payloadKey += ":normal-crease=";
	    payloadKey += std::to_string(normalCreaseAngle);

	    if (payloadKey == presentation.payloadKey) {
		desiredPart = presentation.activePart;
		desiredChannels = presentation.channels;
		desiredGeometryValid = true;
	    } else {
		Obol::PartGeometry geometry;
		const bool progressive =
		    payload->progressiveMesh &&
		    payload->progressiveMesh->isValid();
		bool reusedProgressivePart = false;
		if (progressive) {
		    std::string partKey("compact-progressive-asset:");
		    partKey.reserve(partKey.size() +
			static_cast<size_t>(payload->cacheKey.getLength()) +
			24u);
		    partKey += payload->cacheKey.getString();
		    partKey += ':';
		    partKey += std::to_string(sourceDrawMode);
		    desiredPart = Obol::CadIdBuilder::hash128(partKey);
		}
		const int wire =
		    sourceDrawMode == BOBOL_LOD_DRAW_WIRE || hiddenLine;
		/* SHADED_BOTS creates source-mesh requests only for BoTs.
		 * Those payloads are shaded just like mode 2 payloads. */
		const int shaded =
		    sourceDrawMode == BOBOL_LOD_DRAW_SHADED_BOTS ||
		    sourceDrawMode == BOBOL_LOD_DRAW_SHADED ||
		    hiddenLine;
		std::shared_ptr<const Obol::PartGeometry> preparedGeometry;
		if (progressive &&
		    normalStyle == BObolViewLodState::NORMAL_AUTHORED &&
		    payload->preparedCadGeometry &&
		    payload->preparedCadGeometryRevision ==
			payload->progressiveMesh->revision() &&
		    cad_progressive_geometry_can_draw(
			payload->preparedCadGeometry.get(), wire, shaded,
			desiredActiveLevel)) {
		    preparedGeometry = payload->preparedCadGeometry;
		}

		/*
		 * A newly resident tail must not force a full CPU/GPU realization
		 * while the view is asking for a coarser prefix.  Keep presenting
		 * the older cumulative part when it already covers the active
		 * level.  Leave presentation.payloadKey at the older resident
		 * revision so a later stable, richer request will install the tail.
		 */
		if (progressive && presentation.activePart == desiredPart &&
		    cad_progressive_geometry_can_draw(
			assembly->partGeometry(desiredPart), wire, shaded,
			desiredActiveLevel)) {
		    desiredGeometryValid = true;
		    reusedProgressivePart = true;
		    desiredChannels =
			assembly->compactPartChannels[desiredPart];
		    payloadKey = presentation.payloadKey;
		} else if (progressive &&
		    progressivePartsUpdated.find(desiredPart) !=
			progressivePartsUpdated.end()) {
		    desiredGeometryValid = true;
		    reusedProgressivePart = true;
		    desiredChannels =
			assembly->compactPartChannels[desiredPart];
		} else if (payload->resultKind == BOBOL_LOD_RESULT_AABB ||
		    payload->resultKind == BOBOL_LOD_RESULT_PROXY) {
		    desiredGeometryValid =
			cad_proxy_part_geometry(payload->proxy, geometry) != 0;
		} else if (preparedGeometry) {
		    desiredGeometryValid = true;
		} else {
		    desiredGeometryValid = progressive ?
			cad_progressive_mesh_part_geometry(
			    payload->progressiveMesh, wire, shaded,
			    payload->shadedCullBackfaces, geometry,
			    normalStyle, normalCreaseAngle) != 0 :
			cad_mesh_payload_part_geometry(payload->mesh,
			    wire, shaded, payload->shadedCullBackfaces,
			    geometry, normalStyle, normalCreaseAngle) != 0;
		}
		if (desiredGeometryValid && !reusedProgressivePart) {
		    if (!progressive) {
			std::string partKey("compact-lod:");
			partKey.reserve(partKey.size() +
			    static_cast<size_t>(occurrenceKey.getLength()) +
			    payloadKey.size() + 2u);
			partKey += occurrenceKey.getString();
			partKey += ':';
			partKey += payloadKey;
			desiredPart = Obol::CadIdBuilder::hash128(partKey);
		    }
		    const Obol::PartGeometry *channelGeometry =
			preparedGeometry ? preparedGeometry.get() : &geometry;
		    desiredChannels =
			(channelGeometry->wire.has_value() ? 1u : 0u) |
			(channelGeometry->shaded.has_value() ? 2u : 0u) |
			(channelGeometry->points.has_value() ? 4u : 0u);
		    if (!progressive ||
			progressivePartsUpdated.insert(desiredPart).second) {
			if (preparedGeometry) {
			    Obol::SharedPartUpdate part;
			    part.part = desiredPart;
			    part.geometry = std::move(preparedGeometry);
			    /*
			     * Every immutable PoP generation carries the same
			     * full quantization-domain bounds.  Let Obol validate
			     * and retain occurrence bounds/BVH state while merely
			     * swapping the renderer geometry handle.
			     */
			    part.preservesBounds = true;
			    lodSharedParts.push_back(std::move(part));
			} else {
			    /*
			     * Publish every compact geometry generation through
			     * the immutable shared route.  The legacy PartUpdate
			     * API copies PartGeometry again when it installs the
			     * record; wrapping this completed private value by
			     * move keeps even proxy and non-progressive fallback
			     * publication zero-copy at the assembly boundary.
			     */
			    Obol::SharedPartUpdate part;
			    part.part = desiredPart;
			    part.geometry =
				std::make_shared<const Obol::PartGeometry>(
				    std::move(geometry));
			    lodSharedParts.push_back(std::move(part));
			}
		    }
		    assembly->compactPartChannels[desiredPart] =
			desiredChannels;
		    assembly->compactLodParts.insert(desiredPart);
		}
	    }
	}
	if (!desiredGeometryValid)
	{
	    /* A rejected or stale LoD payload must fall back to the authored
	     * structural part as one coherent record.  Keeping the hashed LoD
	     * part id after geometry conversion failed leaves the instance
	     * referencing a part that was never inserted and makes it vanish. */
	    payloadKey.clear();
	    desiredPart = entry.part;
	    desiredChannels = assembly->compactPartChannels[entry.part];
	    desiredActiveLevel = -1;
	}

	/* Structural proxies and source meshes generally do not share a local
	 * coordinate system.  Cold coverage may use one normalized unit box
	 * transformed by geometryTransform, while the cached PoP arrays remain
	 * in the BoT asset's original coordinates.  Swapping only the part id
	 * therefore draws a valid mesh with the proxy transform—often far
	 * outside the view or collapsed to an apparently empty wireframe.
	 * Choose the transform with the chosen representation and update it in
	 * the same instance mutation as the part. */
	const SbMatrix desiredLocalToRoot =
	    desiredActiveLevel >= 0 && entry.sourceMeshRequestValid ?
		compact_mesh_asset_matrix(this, entry) : entry.localToSource;

	const bool partChanged = presentation.activePart != desiredPart ||
	    presentation.payloadKey != payloadKey;
	const bool levelChanged =
	    presentation.activeLevel != desiredActiveLevel;
	const bool appearanceChanged = !reset &&
	    presentation.appearanceRevision != entry.appearanceRevision;
	const bool selectionChanged = !reset &&
	    presentation.selectionRevision != entry.selectionRevision;
	const bool recordChanged = !reset &&
	    (partChanged ||
	     presentation.geometryRevision != entry.geometryRevision ||
	     presentation.placementRevision != entry.placementRevision);
	if (partChanged) {
	    if (presentation.activePart != entry.part &&
		presentation.activePart != desiredPart)
		partsToRemove.insert(presentation.activePart);
	    presentation.payloadKey = payloadKey;
	    presentation.activePart = desiredPart;
	    presentation.channels = desiredChannels;
	}
	if (!presentationWasPresent || partChanged) {
	    if (presentationWasPresent) {
		const auto previousReference =
		    assembly->compactActivePartReferences.find(previousPart);
		if (previousReference !=
			assembly->compactActivePartReferences.end()) {
		    if (previousReference->second > 1)
			previousReference->second--;
		    else
			assembly->compactActivePartReferences.erase(
			    previousReference);
		}
	    }
	    assembly->compactActivePartReferences[desiredPart]++;
	}
	if (reset && i < instances.size()) {
	    instances[i].record.part = desiredPart;
	    instances[i].record.localToRoot = desiredLocalToRoot;
	    instances[i].record.lodLevel = desiredActiveLevel >= 0 ?
		static_cast<uint8_t>(std::min(15, desiredActiveLevel)) : 255;
	} else if (recordChanged &&
	    i < this->d->compactIndex->instances.size()) {
	    Obol::InstanceUpdate update =
		this->d->compactIndex->instances[i];
	    update.record.part = desiredPart;
	    update.record.localToRoot = desiredLocalToRoot;
	    update.record.lodLevel = desiredActiveLevel >= 0 ?
		static_cast<uint8_t>(std::min(15, desiredActiveLevel)) : 255;
	    instances.push_back(std::move(update));
	} else if (appearanceChanged || selectionChanged) {
	    Obol::InstanceStyleUpdate update;
	    update.instance = entry.instance;
	    update.style = entry.style;
	    instanceStyles.push_back(update);
	} else if (levelChanged) {
	    Obol::InstanceLodUpdate update;
	    update.instance = entry.instance;
	    update.lodLevel = desiredActiveLevel >= 0 ?
		static_cast<uint8_t>(std::min(15, desiredActiveLevel)) : 255;
	    lodLevelUpdates.push_back(update);
	}
	presentation.activeLevel = desiredActiveLevel;

	if (incrementalUpdate) {
	    const bool previousWire = previousChannels & (1u | 4u);
	    const bool currentWire = presentation.channels & (1u | 4u);
	    const bool previousShaded = previousChannels & 2u;
	    const bool currentShaded = presentation.channels & 2u;
	    if (previousWire != currentWire)
		wireCount = currentWire ? wireCount + 1 :
		    (wireCount > 0 ? wireCount - 1 : 0);
	    if (previousShaded != currentShaded)
		shadedCount = currentShaded ? shadedCount + 1 :
		    (shadedCount > 0 ? shadedCount - 1 : 0);
	} else {
	    if (presentation.channels & (1u | 4u))
		wireCount++;
	    if (presentation.channels & 2u)
		shadedCount++;
	}

	if (!reset &&
	    (presentation.visibilityRevision != entry.visibilityRevision ||
	     presentation.selectionRevision != entry.selectionRevision)) {
	    instanceSetsChanged = true;
	}
	presentation.geometryRevision = entry.geometryRevision;
	presentation.appearanceRevision = entry.appearanceRevision;
	presentation.placementRevision = entry.placementRevision;
	presentation.visibilityRevision = entry.visibilityRevision;
	presentation.selectionRevision = entry.selectionRevision;
    }
    const int64_t presentationPlanned =
	presentationTiming ? bu_gettime() : 0;

    /* Progressive assets are shared by every occurrence of the same source
     * geometry.  One occurrence changing back to its base part must not
     * remove a retained part that another occurrence still references. */
    if (!partsToRemove.empty()) {
	for (auto part = partsToRemove.begin();
	     part != partsToRemove.end();) {
	    const auto references =
		assembly->compactActivePartReferences.find(*part);
	    if (references !=
		    assembly->compactActivePartReferences.end() &&
		references->second > 0)
		part = partsToRemove.erase(part);
	    else
		++part;
	}
    }

    int presentationDrawMode = SoCADAssembly::WIREFRAME;
    if (hiddenLine)
	presentationDrawMode = SoCADAssembly::HIDDEN_LINE;
    else if (shadedCount > 0 && wireCount > 0)
	presentationDrawMode = SoCADAssembly::SHADED_WITH_EDGES;
    else if (shadedCount > 0)
	presentationDrawMode = SoCADAssembly::SHADED;
    const bool drawModeChanged =
	assembly->drawMode.getValue() != presentationDrawMode;

    const bool structuralChanged = reset || !lodSharedParts.empty() ||
	!baseParts.empty() || !instances.empty() || !instanceStyles.empty() ||
	!partsToRemove.empty() ||
	instanceSetsChanged || drawModeChanged;
    if (structuralChanged) {
	assembly->beginUpdate();
	if (reset) {
	    assembly->clear();
	    assembly->clearSemanticMap();
	    for (const BObolCompactInstanceEntry &entry :
		 this->d->compactIndex->entries) {
		SoBRLCadAssembly::InstanceSemantic semantic = entry.semantic;
		semantic.sourceInstanceKey = compact_instance_identity(entry);
		assembly->setInstanceSemantic(entry.instance, semantic);
	    }
	}
	assembly->upsertSharedParts(baseParts);
	assembly->upsertSharedParts(lodSharedParts);
	assembly->upsertInstances(instances);
	assembly->updateInstanceStyles(instanceStyles);
	std::vector<Obol::PartId> removedParts;
	removedParts.reserve(partsToRemove.size());
	for (const Obol::PartId &part : partsToRemove) {
	    removedParts.push_back(part);
	    assembly->compactPartChannels.erase(part);
	    assembly->compactLodParts.erase(part);
	}
	if (!reset) {
	    for (const size_t i : structuralEntryIndices) {
		if (i >= this->d->compactIndex->entries.size())
		    continue;
		const BObolCompactInstanceEntry &entry =
		    this->d->compactIndex->entries[i];
		SoBRLCadAssembly::InstanceSemantic semantic = entry.semantic;
		semantic.sourceInstanceKey = compact_instance_identity(entry);
		assembly->setInstanceSemantic(entry.instance, semantic);
	    }
	}
	assembly->removeParts(removedParts);
	if (instanceSetsChanged) {
	    assembly->setHiddenInstances(
		this->d->compactIndex->hiddenInstances);
	    assembly->setSelectedInstances(
		this->d->compactIndex->selectedInstances);
	    assembly->setUnpickableInstances(
		this->d->compactIndex->unpickableInstances);
	}
	assembly->drawMode = presentationDrawMode;
	assembly->pickMode = SoCADAssembly::PICK_AUTO;
	assembly->endUpdate();
    }
    if (!lodLevelUpdates.empty())
	assembly->updateInstanceLodLevels(lodLevelUpdates);

    /*
     * Compact presentation is intentionally sparse, so an intended active
     * part and the retained Obol instance must never drift silently.  Keep
     * this full invariant audit behind an opt-in diagnostic: it is useful for
     * validating streaming box-to-mesh waves, but an O(N) scan does not
     * belong in the normal publication path for very large scenes.
     */
    const char *presentationDebug =
	getenv("BOBOL_COMPACT_PRESENTATION_DEBUG");
    if (presentationDebug && presentationDebug[0] &&
	!BU_STR_EQUAL(presentationDebug, "0")) {
	size_t missingInstances = 0;
	size_t mismatchedParts = 0;
	size_t retainedStructural = 0;
	size_t structuralLodBacked = 0;
	size_t structuralSourceRequests = 0;
	size_t structuralEntryGeometry = 0;
	size_t structuralGeometryPointerMismatch = 0;
	std::map<std::string, size_t> structuralKinds;
	std::map<std::string, size_t> structuralTypes;
	for (const auto &item : assembly->compactInstancePresentations) {
	    const std::optional<Obol::InstanceRecord> record =
		assembly->getInstanceRecord(item.first);
	    if (!record) {
		missingInstances++;
		continue;
	    }
	    if (!(record->part == item.second.activePart))
		mismatchedParts++;
	    const Obol::PartGeometry *geometry =
		assembly->partGeometry(record->part);
	    if (geometry && geometry->structuralProxy) {
		retainedStructural++;
		const auto entryIndex =
		    this->d->compactIndex->entryIndex.find(item.first);
		if (entryIndex !=
			this->d->compactIndex->entryIndex.end() &&
		    entryIndex->second <
			this->d->compactIndex->entries.size()) {
		    const BObolCompactInstanceEntry &entry =
			this->d->compactIndex->entries[entryIndex->second];
		    if (entry.geometry && entry.geometry->structuralProxy)
			structuralEntryGeometry++;
		    if (!entry.geometry ||
			entry.geometry.get() != geometry)
			structuralGeometryPointerMismatch++;
		    if (entry.lodBacked)
			structuralLodBacked++;
		    if (entry.sourceMeshRequestValid)
			structuralSourceRequests++;
		    structuralKinds[
			entry.shapeSummary.geometryKind.getString()]++;
		    structuralTypes[
			entry.shapeSummary.sourceType.getString()]++;
		}
	    }
	}
	bu_log("BObol compact presentation invariant presentations=%zu "
	       "retained_instances=%zu retained_structural=%zu "
	       "structural_lod=%zu structural_requests=%zu "
	       "structural_entry_geometry=%zu structural_pointer_mismatch=%zu "
	       "missing_instances=%zu mismatched_parts=%zu "
	       "updates=%zu lod_updates=%zu reset=%d\n",
	       assembly->compactInstancePresentations.size(),
	       assembly->instanceCount(), retainedStructural,
	       structuralLodBacked, structuralSourceRequests,
	       structuralEntryGeometry, structuralGeometryPointerMismatch,
	       missingInstances, mismatchedParts, instances.size(),
	       lodLevelUpdates.size(), reset ? 1 : 0);
	if (retainedStructural) {
	    for (const auto &kind : structuralKinds)
		bu_log("  structural kind=%s count=%zu\n",
		    kind.first.c_str(), kind.second);
	    for (const auto &type : structuralTypes)
		bu_log("  structural type=%s count=%zu\n",
		    type.first.c_str(), type.second);
	}
    }

    assembly->compactPresentationInitialized = TRUE;
    assembly->compactPresentationIndex = this->d->compactIndex;
    assembly->compactPresentationSourceRevision =
	this->sourceRevision.getValue();
    assembly->compactPresentationInputsRevision =
	this->inputsRevision.getValue();
    assembly->compactPresentationCadBatchRevision =
	this->cadBatchRevisionGet();
    assembly->compactPresentationPayloadRevision = payloadRevision;
    assembly->compactPresentationDrawMode = sourceDrawMode;
    assembly->compactWirePresentationCount = wireCount;
    assembly->compactShadedPresentationCount = shadedCount;
    viewState->acknowledgeCadOccurrenceChanges(this, payloadRevision);
    if (presentationTiming) {
	const int64_t presentationCompleted = bu_gettime();
	const int64_t total =
	    presentationCompleted - presentationStarted;
	if (total >= 10000) {
	    bu_log("BObol compact presentation total=%" PRId64 "us "
		   "index=%" PRId64 "us plan=%" PRId64 "us "
		   "publish=%" PRId64 "us reset=%d full_resync=%d "
		   "sparse_payload=%d cad_batch_changed=%d "
		   "changed_keys=%zu structural_entries=%zu entries=%zu "
		   "payload_input=%zu shared_parts=%zu "
		   "instance_updates=%zu lod_updates=%zu\n",
		   total,
		   presentationIndexed - presentationStarted,
		   presentationPlanned - presentationIndexed,
		   presentationCompleted - presentationPlanned,
		   reset ? 1 : 0, fullPayloadResync ? 1 : 0,
		   sparsePayloadUpdate ? 1 : 0,
		   cadBatchChanged ? 1 : 0,
		   changedOccurrenceKeys.size(),
		   structuralEntryIndices.size(), entryIndices.size(),
		   payloadInput->size(), lodSharedParts.size(),
		   instances.size(),
		   lodLevelUpdates.size());
	}
    }
    return assembly;
}

static SoBRLCadAssembly *
cad_view_lod_assembly_for_action(SoAction *action,
				 const SoBRLDatabaseSource *source)
{
    const BObolViewLodState *viewState =
	bobol_view_lod_state_for_action(action);
    if (source && source->isCompactOccurrenceRegistry()) {
	if (viewState) {
	    SoBRLCadAssembly *current =
		source->currentCompactViewLodAssembly(viewState);
	    if (current)
		return current;
	    std::vector<const BObolViewLodState::CadPayload *> payloads;
	    return source->compactViewLodAssembly(payloads, viewState);
	}
	/* Aggregate sources require an exact source binding.  The legacy
	 * path/name lookup may otherwise attach a sibling source's payload. */
	return NULL;
    }
    const BObolViewLodState::CadPayload *payload =
	bobol_view_lod_cad_for_action(action, source);
    return cad_view_lod_assembly(source, payload, viewState);
}

static void
cad_add_vlist_instance(cad_build_data &data,
		       SoBRLVListShape *shape,
		       const SbMatrix &localMatrix)
{
    if (!shape)
	return;
    Obol::PartGeometry geometry;
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
    Obol::InstanceId instanceId = Obol::CadIdBuilder::hash128(instanceKey);
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
    if (data.storeSemantics)
	data.assembly->setInstanceSemantic(instanceId, cad_vlist_semantic(shape));
    data.wireCount++;
}

static void
cad_add_mesh_instance(cad_build_data &data,
		      SoBRLMeshShape *shape,
		      const SbMatrix &localMatrix)
{
    if (!shape)
	return;
    Obol::PartGeometry geometry;
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
    Obol::InstanceId instanceId = Obol::CadIdBuilder::hash128(instanceKey);
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
    if (data.storeSemantics)
	data.assembly->setInstanceSemantic(instanceId, cad_mesh_semantic(shape));
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

static void
compact_add_part_if_needed(BObolCompactInstanceIndex &index,
			   const std::string &partKey,
			   const std::shared_ptr<const Obol::PartGeometry> &geometry,
			   Obol::PartId &partId)
{
    std::map<std::string, Obol::PartId>::iterator found =
	index.partIdByKey.find(partKey);
    if (found != index.partIdByKey.end()) {
	partId = found->second;
	return;
    }

    partId = Obol::CadIdBuilder::hash128(partKey);
    index.partIdByKey[partKey] = partId;
    BObolCompactPartReference part;
    part.part = partId;
    part.geometry = geometry;
    index.parts.push_back(part);
}

static int
compact_add_geometry_part_if_needed(BObolCompactInstanceIndex &index,
	const char *kind,
	const std::shared_ptr<const Obol::PartGeometry> &geometry,
	Obol::PartId &partId)
{
    if (!geometry)
	return 0;
    const auto found = index.partIdByGeometry.find(geometry.get());
    if (found != index.partIdByGeometry.end()) {
	partId = found->second;
	return 1;
    }

    std::string partKey;
    if (!cad_part_key_for_geometry(kind, *geometry, partKey))
	return 0;
    compact_add_part_if_needed(index, partKey, geometry, partId);
    index.partIdByGeometry[geometry.get()] = partId;
    return 1;
}

static void compact_sync_entry_from_source(BObolCompactInstanceEntry &entry,
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
	int &unsupported)
{
    const std::shared_ptr<const Obol::PartGeometry> &geometry =
	input.occurrence.geometry;
    if (!source || !geometry) {
	unsupported = 1;
	return;
    }

    const char *partKind = geometry->shaded ? "mesh" :
	(geometry->points && !geometry->wire ? "point" : "wire");
    Obol::PartId partId;
    if (!compact_add_geometry_part_if_needed(index,
	    partKind, geometry, partId)) {
	unsupported = 1;
	return;
    }

    const char *path = input.occurrence.summary.path.getString();
    if (!path || !path[0])
	path = source->path.getValue().getString();
    SbMatrix geometryToSource = input.occurrence.geometryTransform;
    geometryToSource.multRight(input.occurrence.localTransform);
    const SbMatrix matrix = cad_instance_matrix(source, geometryToSource);
    const std::string instanceKey = cad_instance_key(source, path, ordinal++);
    const Obol::InstanceId instanceId =
	Obol::CadIdBuilder::hash128(instanceKey);

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
    entry.lodBacked = input.occurrence.lodBacked;
    entry.sourceMeshRequestValid = input.occurrence.sourceMeshRequestValid;
    if (entry.sourceMeshRequestValid) {
	entry.sourceMeshRequest = input.occurrence.sourceMeshRequest;
	index.sourceMeshRequestCount++;
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
    entry.highlighted = input.occurrence.summary.highlighted;
    entry.shapeSummary = input.occurrence.summary;
    entry.occurrenceIndex = input.occurrence.occurrenceIndex;
    entry.booleanOperation = input.occurrence.booleanOperation;
    compact_sync_entry_from_source(entry, source);
    if (input.stylesValid) {
	entry.normalStyle = input.normalStyle;
	entry.selectedStyle = input.selectedStyle;
	entry.highlightedStyle = input.highlightedStyle;
    }
    if (input.dashed) {
	entry.normalStyle.linePattern = 0xcf33u;
	entry.selectedStyle.linePattern = 0xcf33u;
	entry.highlightedStyle.linePattern = 0xcf33u;
    }
    entry.style = compact_effective_style(entry);
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
    record.style = entry.style;

    Obol::InstanceUpdate update;
    update.instance = instanceId;
    update.record = record;
    index.instances.push_back(update);
    if (!entry.visible)
	index.hiddenInstances.push_back(instanceId);
    if (entry.selected)
	index.selectedInstances.push_back(instanceId);
    if (!entry.selectable)
	index.unpickableInstances.push_back(instanceId);
    index.entries.push_back(entry);
    if (entry.meshGeometry)
	index.shadedCount++;
    if (entry.wireGeometry || entry.pointGeometry)
	index.wireCount++;
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

static void
compact_index_bounds_remove(BObolCompactInstanceIndex &index,
	BObolCompactInstanceEntry &entry)
{
    if (entry.sourceBounds.isEmpty())
	return;

    if (!index.sourceBoundsDirty && !index.sourceBounds.isEmpty()) {
	const SbVec3f entryMinimum = entry.sourceBounds.getMin();
	const SbVec3f entryMaximum = entry.sourceBounds.getMax();
	const SbVec3f aggregateMinimum = index.sourceBounds.getMin();
	const SbVec3f aggregateMaximum = index.sourceBounds.getMax();
	for (int axis = 0; axis < 3; axis++) {
	    if (entryMinimum[axis] <= aggregateMinimum[axis] ||
		entryMaximum[axis] >= aggregateMaximum[axis]) {
		index.sourceBoundsDirty = true;
		break;
	    }
	}
    }
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

static void
compact_merge_runtime_state(const BObolCompactInstanceIndex *current,
    BObolCompactInstanceIndex *next)
{
    if (!next)
	return;
    if (!current) {
	compact_rebuild_entry_index(*next);
	return;
    }
    std::unordered_map<Obol::InstanceId,
	const BObolCompactInstanceEntry *, std::hash<Obol::InstanceId>> old;
    std::unordered_map<std::string, const BObolCompactInstanceEntry *>
	oldByPath;
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

    next->hiddenInstances.clear();
    next->selectedInstances.clear();
    next->unpickableInstances.clear();
    std::unordered_set<Obol::InstanceId, std::hash<Obol::InstanceId>>
	assignedInstances;
    assignedInstances.reserve(next->entries.size());
    for (size_t i = 0; i < next->entries.size(); i++) {
	BObolCompactInstanceEntry &entry = next->entries[i];
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
		if (i < next->instances.size()) {
		    next->instances[i].instance = previous.instance;
		}
	    }
	    entry.authoredVisible = previous.authoredVisible;
	    entry.visible = previous.visible;
	    entry.selectable = previous.selectable;
	    entry.selected = previous.selected;
	    entry.highlighted = previous.highlighted;
	    entry.geometryRevision = previous.geometryRevision;
	    if (entry.part != previous.part)
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
	}
	assignedInstances.insert(entry.instance);
	entry.style = compact_effective_style(entry);
	if (previousEntry && !compact_appearance_equal(entry, *previousEntry))
	    entry.appearanceRevision = compact_next_revision(
		entry.appearanceRevision);
	if (i < next->instances.size())
	    next->instances[i].record.style = entry.style;
	compact_sync_shape_summary(entry);
	if (!entry.visible)
	    next->hiddenInstances.push_back(entry.instance);
	if (entry.selected)
	    next->selectedInstances.push_back(entry.instance);
	if (!entry.selectable)
	    next->unpickableInstances.push_back(entry.instance);
    }
    compact_rebuild_entry_index(*next);
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

    BObolCompactInstanceIndex *next = new BObolCompactInstanceIndex;
    compact_occurrence_build input;
    input.occurrence = occurrence;
    input.semantic = compact_semantic_from_summary(occurrence.summary);
    int ordinal = 0;
    int unsupported = 0;
    compact_add_occurrence(this, *next, input, ordinal, unsupported);
    if (unsupported || next->entries.empty()) {
	delete next;
	return 0;
    }

    this->installCompactInstanceIndex(next, FALSE);
    remove_non_auxiliary_children(this);
    this->markCompiledAssemblyDirty();

    SbMatrix geometryToSource = occurrence.geometryTransform;
    geometryToSource.multRight(occurrence.localTransform);
    const SbBox3f bounds = database_source_transform_bounds(
	compact_part_geometry_bounds(occurrence.geometry), geometryToSource);
    if (!bounds.isEmpty())
	(void)this->setSourceBoundsState(TRUE, bounds.getMin(), bounds.getMax());
    else
	this->clearSourceBounds();

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

    BObolCompactInstanceIndex *next = new BObolCompactInstanceIndex;
    int ordinal = 0;
    int unsupported = 0;
    SbBox3f aggregateBounds;
    aggregateBounds.makeEmpty();
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
	SbMatrix geometryToSource = occurrence.geometryTransform;
	geometryToSource.multRight(occurrence.localTransform);
	const SbBox3f bounds = database_source_transform_bounds(
	    compact_part_geometry_bounds(occurrence.geometry),
	    geometryToSource);
	if (!bounds.isEmpty())
	    aggregateBounds.extendBy(bounds);
    }
    if (unsupported || next->entries.size() != occurrences.size()) {
	delete next;
	return 0;
    }

    this->installCompactInstanceIndex(next, TRUE);
    remove_non_auxiliary_children(this);
    this->markCompiledAssemblyDirty();
    if (!aggregateBounds.isEmpty())
	(void)this->setSourceBoundsState(TRUE, aggregateBounds.getMin(),
	    aggregateBounds.getMax());
    else
	this->clearSourceBounds();

    bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_SOURCES, 1);
    bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_INSTANCES,
	static_cast<uint64_t>(this->d->compactIndex->entries.size()));
    return this->d->compactIndex->entries.size() >
	static_cast<size_t>(INT_MAX) ? INT_MAX :
	static_cast<int>(this->d->compactIndex->entries.size());
}

/* Populate the derived lookup maps and aggregate bounds for the appended tail
 * [firstNew, entries.size()).  compact_add_occurrence intentionally leaves these
 * to a rebuild at install time; mergeCompactOccurrences must fill them for the
 * new entries only (mirrors compact_rebuild_entry_index's per-entry work) so the
 * live registry stays consistent without an O(N) full rebuild each batch. */
static void
compact_index_append_derived(BObolCompactInstanceIndex &index, size_t firstNew)
{
    for (size_t i = firstNew; i < index.entries.size(); i++) {
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
	const char *sourceName =
	    index.entries[i].semantic.sourceName.getString();
	if (sourceName && sourceName[0])
	    index.entryIndicesBySourceName[sourceName].push_back(i);
	index.partReferenceCounts[index.entries[i].part]++;
	compact_index_bounds_add(index, index.entries[i]);
    }
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

/* Replace one leaf occurrence's drawing data in place: swap its part geometry,
 * transform, flags, and summary for the incoming higher-tier occurrence while
 * preserving the entry's instance identity, semantic key, and runtime
 * (visible/selected/pickable) state.  instances[i] runs parallel to entries[i]
 * (both are appended together by compact_add_occurrence and never reordered), so
 * the instance record is updated at the same index. */
static bool
compact_index_replace_entry_geometry(SoBRLDatabaseSource *source,
	BObolCompactInstanceIndex &index, size_t entryIdx,
	const BObolCompactOccurrence &occurrence)
{
    const std::shared_ptr<const Obol::PartGeometry> &geometry =
	occurrence.geometry;
    if (!source || !geometry || entryIdx >= index.entries.size() ||
	entryIdx >= index.instances.size())
	return false;

    const char *partKind = geometry->shaded ? "mesh" :
	(geometry->points && !geometry->wire ? "point" : "wire");
    Obol::PartId newPartId;
    if (!compact_add_geometry_part_if_needed(index, partKind, geometry,
	    newPartId))
	return false;

    BObolCompactInstanceEntry &entry = index.entries[entryIdx];
    Obol::InstanceUpdate &update = index.instances[entryIdx];

    /* Drop the outgoing entry's bounds/part-ref/kind-count contribution before
     * mutating it, then re-add for the replacement. */
    compact_index_bounds_remove(index, entry);
    std::unordered_map<Obol::PartId, size_t, std::hash<Obol::PartId>>::iterator
	oldRef = index.partReferenceCounts.find(entry.part);
    if (oldRef != index.partReferenceCounts.end() && oldRef->second > 0)
	oldRef->second--;
    if (entry.meshGeometry) {
	if (index.shadedCount > 0)
	    index.shadedCount--;
    }
    if ((entry.wireGeometry || entry.pointGeometry) &&
	index.wireCount > 0) {
	index.wireCount--;
    }

    SbMatrix geometryToSource = occurrence.geometryTransform;
    geometryToSource.multRight(occurrence.localTransform);
    const SbMatrix matrix = cad_instance_matrix(source, geometryToSource);

    entry.part = newPartId;
    entry.geometry = geometry;
    entry.wireGeometry = geometry->wire ? TRUE : FALSE;
    entry.pointGeometry = geometry->points ? TRUE : FALSE;
    entry.meshGeometry = geometry->shaded ? TRUE : FALSE;
    entry.lodBacked = occurrence.lodBacked;
    const SbBool oldSourceMeshRequestValid = entry.sourceMeshRequestValid;
    entry.sourceMeshRequestValid = occurrence.sourceMeshRequestValid;
    if (entry.sourceMeshRequestValid)
	entry.sourceMeshRequest = occurrence.sourceMeshRequest;
    else
	entry.sourceMeshRequest.clear();
    if (oldSourceMeshRequestValid != entry.sourceMeshRequestValid) {
	if (entry.sourceMeshRequestValid)
	    index.sourceMeshRequestCount++;
	else if (index.sourceMeshRequestCount > 0)
	    index.sourceMeshRequestCount--;
    }
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

    update.record.part = newPartId;
    update.record.localToRoot = matrix;
    update.record.style = entry.style;

    index.partReferenceCounts[newPartId]++;
    compact_index_bounds_add(index, entry);
    if (entry.meshGeometry)
	index.shadedCount++;
    if (entry.wireGeometry || entry.pointGeometry)
	index.wireCount++;
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

    BObolCompactInstanceEntry &entry = index.entries[entryIdx];
    bool changed = false;
    if (occurrence.lodBacked && !entry.lodBacked) {
	entry.lodBacked = TRUE;
	changed = true;
    }
    if (occurrence.sourceMeshRequestValid) {
	if (!entry.sourceMeshRequestValid) {
	    entry.sourceMeshRequestValid = TRUE;
	    index.sourceMeshRequestCount++;
	}
	entry.sourceMeshRequest = occurrence.sourceMeshRequest;
	changed = true;
    }
    if (!changed)
	return false;

    /* Source semantics become authoritative with the request, but the
     * currently drawn geometry, proxy transform, placement, part, and runtime
     * visibility/selection state remain untouched.
     *
     * This is a data-management contract change, not a drawable-geometry
     * change.  Do not advance geometryRevision: doing so queues a redundant
     * full InstanceRecord update for every enriched leaf even though the
     * retained part, transform, bounds, style, and selection state are
     * byte-for-byte unchanged.  markDisplayMeshLodDirty below independently
     * makes the submit planner observe the new source request. */
    entry.shapeSummary.shapeKind = occurrence.summary.shapeKind;
    entry.shapeSummary.sourceType = occurrence.summary.sourceType;
    entry.shapeSummary.sourceId = occurrence.summary.sourceId;
    compact_sync_shape_summary(entry);
    return true;
}

int
SoBRLDatabaseSource::mergeCompactOccurrences(
    const std::vector<BObolCompactOccurrence> &occurrences)
{
    if (occurrences.empty())
	return 0;

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
    std::vector<size_t> changedEntries;
    changedEntries.reserve(occurrences.size());
    for (const BObolCompactOccurrence &occurrence : occurrences) {
	if (!occurrence.geometry)
	    continue;
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
		 evolvingOverview);
	    if (newTier < oldTier ||
		(newTier == oldTier && !richerDataContract))
		continue;
	    const bool merged = newTier == oldTier && !evolvingOverview ?
		compact_index_merge_source_contract(index, found->second,
		    occurrence) :
		compact_index_replace_entry_geometry(this, index,
		    found->second, occurrence);
	    if (merged) {
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
	compact_add_occurrence(this, index, input, ordinal, unsupported);
	if (unsupported || index.entries.size() == before)
	    continue;
	if (p && p[0])
	    index.entryIndexByPath[p] = index.entries.size() - 1;
	changed++;
	changedEntries.push_back(index.entries.size() - 1);
	if (occurrence.sourceMeshRequestValid)
	    mergedSourceContract = true;
    }
    if (!changed)
	return 0;

    if (index.entries.size() > firstNew)
	compact_index_append_derived(index, firstNew);
    if (!this->d->compactExpectedInstanceCount)
	this->d->compactExpectedInstanceCount = index.entries.size();
    /*
     * Newly appended occurrences inherit authored visibility and start
     * unselected.  Rebuilding N-entry masks after every streamed batch when
     * neither a visibility frontier nor selected paths exist revisits the
     * complete accumulated scene and turns otherwise append-only realization
     * into O(N^2) GUI-thread work.
     *
     * Active frontiers still need to classify newly arriving paths.  Their
     * existing reapply routines remain authoritative until they acquire a
     * sparse changed-entry form.
     */
    if (this->d->compactVisibilityFrontierActive)
	(void)this->reapplyCompactInstanceVisibilityFrontier(firstNew);
    if (!this->d->compactSelectedPaths.empty())
	(void)this->reapplyCompactInstanceSelectedPaths(firstNew);
    this->d->compactOccurrenceRegistry = TRUE;
    this->markCompiledAssemblyDirty();
    /*
     * Streaming is append/upgrade-only.  Preserve that exact sparse
     * mutation contract for both the retained assembly and view-LoD
     * consumers instead of invalidating their complete source cursors.
     * The former source-wide dirty marks made a one-leaf append force every
     * view to revisit all previously published leaves, which is quadratic
     * during cold population and can starve the unvisited tail.
     */
    std::sort(changedEntries.begin(), changedEntries.end());
    changedEntries.erase(std::unique(changedEntries.begin(),
	changedEntries.end()), changedEntries.end());
    this->markCadBatchDirty(changedEntries);
    this->markDisplayMeshLodDirty(changedEntries);
    if (mergedSourceContract) {
	this->d->displayMeshLodContractRevisionValid = TRUE;
	this->d->displayMeshLodContractSourceRevision =
	    this->sourceRevision.getValue();
	this->d->displayMeshLodContractInputsRevision =
	    this->inputsRevision.getValue();
    }

    SbBox3f bounds;
    if (compact_index_source_bounds(index, bounds) && !bounds.isEmpty())
	(void)this->setSourceBoundsState(TRUE, bounds.getMin(),
	    bounds.getMax());
    return changed;
}

void
SoBRLDatabaseSource::reserveCompactOccurrenceCapacity(size_t expectedCount)
{
    if (!expectedCount)
	return;
    this->d->compactExpectedInstanceCount = std::max(
	this->d->compactExpectedInstanceCount, expectedCount);
    if (!this->d->compactIndex) {
	this->d->compactIndex = new BObolCompactInstanceIndex;
	this->d->compactIndexActive = TRUE;
	this->d->compactOccurrenceRegistry = TRUE;
    }
    BObolCompactInstanceIndex &index = *this->d->compactIndex;
    const size_t capacity = expectedCount < SIZE_MAX ?
	expectedCount + 1 : expectedCount;
    if (capacity <= index.entries.capacity())
	return;

    index.entries.reserve(capacity);
    index.instances.reserve(capacity);
    index.parts.reserve(capacity);
    /*
     * Do not eagerly allocate every hash table to the full manifest count.
     * Ten 50k bucket arrays produced a 70+ ms UI-thread pause before the first
     * batch.  Reserving the relocation-sensitive contiguous records prevents
     * the observed late-stream copy spike; the associative indices then grow
     * with the bounded publication batches.
     */
}


int
SoBRLDatabaseSource::adoptCompactOccurrencesFrom(
    const SoBRLDatabaseSource *source)
{
    if (!source || source == this || !source->d->compactIndex ||
	source->d->compactIndex->entries.empty())
	return 0;

    SbMatrix placementInverse = SbMatrix::identity();
    if (this->drawMatrixValid.getValue())
	placementInverse = this->drawMatrix.getValue().inverse();

    std::vector<BObolCompactOccurrence> occurrences;
    occurrences.reserve(source->d->compactIndex->entries.size());
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
	occurrence.lodBacked = entry.lodBacked;
	occurrence.sourceMeshRequestValid = entry.sourceMeshRequestValid;
	occurrence.sourceMeshRequest = entry.sourceMeshRequest;
	occurrence.occurrenceIndex = entry.occurrenceIndex;
	occurrence.booleanOperation = entry.booleanOperation;
	occurrences.push_back(occurrence);
    }
    return this->mergeCompactOccurrences(occurrences);
}

static bool
compact_retarget_path_component(std::string &path,
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


int
SoBRLDatabaseSource::retargetCompactOccurrencePaths(
    const char *oldPath, const char *newPath)
{
    if (!this->d->compactIndex || !oldPath || !oldPath[0] ||
	!newPath || !newPath[0])
	return 0;

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
    for (size_t i = 0; i < this->d->compactIndex->entries.size(); i++) {
	BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[i];
	std::string semanticPath = entry.semantic.path.getString();
	if (!compact_retarget_path_component(semanticPath, oldComponent,
		newComponent))
	    continue;

	entry.semantic.path = semanticPath.c_str();
	if (database_source_db_path_without_instance_suffixes(
		entry.semantic.sourceName.getString()) == oldComponent)
	    entry.semantic.sourceName = newComponent.c_str();
	if (entry.sourceMeshRequestValid) {
	    std::string requestPath =
		entry.sourceMeshRequest.path.getString();
	    (void)compact_retarget_path_component(requestPath, oldComponent,
		newComponent);
	    entry.sourceMeshRequest.path = requestPath.c_str();
	    std::string assetPath =
		entry.sourceMeshRequest.meshAssetPath.getString();
	    (void)compact_retarget_path_component(assetPath, oldComponent,
		newComponent);
	    entry.sourceMeshRequest.meshAssetPath = assetPath.c_str();
	    if (BU_STR_EQUAL(
		    entry.sourceMeshRequest.sourceName.getString(),
		    oldComponent.c_str()))
		entry.sourceMeshRequest.sourceName = newComponent.c_str();
	    if (BU_STR_EQUAL(
		    entry.sourceMeshRequest.meshAssetName.getString(),
		    oldComponent.c_str()))
		entry.sourceMeshRequest.meshAssetName = newComponent.c_str();
	}
	if (i < this->d->compactIndex->instances.size())
	    this->d->compactIndex->instances[i].record.childName = semanticPath;
	compact_sync_shape_summary_state(entry);
	changed++;
    }
    if (!changed)
	return 0;

    compact_rebuild_entry_index(*this->d->compactIndex);
    this->markCompiledAssemblyDirty();
    this->markCadBatchDirty();
    this->markDisplayMeshLodDirty();
    this->touch();
    return changed;
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

static uint64_t
compact_state_signature(const std::vector<Obol::InstanceId> &instances)
{
    uint64_t signature = 1469598103934665603ULL;
    for (const Obol::InstanceId &instance : instances) {
	signature ^= instance.w0;
	signature *= 1099511628211ULL;
	signature ^= instance.w1;
	signature *= 1099511628211ULL;
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
    int drawMode = SoCADAssembly::WIREFRAME;
    if (source_record_draw_mode(source) == BOBOL_LOD_DRAW_HIDDEN_LINE)
	drawMode = SoCADAssembly::HIDDEN_LINE;
    else if (index->shadedCount > 0 && index->wireCount > 0)
	drawMode = SoCADAssembly::SHADED_WITH_EDGES;
    else if (index->shadedCount > 0)
	drawMode = SoCADAssembly::SHADED;
    if (assembly->drawMode.getValue() != drawMode)
	assembly->drawMode = drawMode;
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

int
SoBRLDatabaseSource::syncCompiledAssembly(void)
{
    const SbUniqueId sourceNodeId = this->getNodeId();
    if (!this->d->compiledAssemblyDirty &&
	this->d->compiledAssemblyNodeId == sourceNodeId)
	return this->d->compiledAssemblyActive ? 1 : 0;
    this->d->compiledAssemblyActive = FALSE;

    const SbBool hasCompactPayload = this->d->compactIndexActive &&
	this->d->compactIndex && !this->d->compactIndex->instances.empty();
    if (!this->visible.getValue() ||
	this->auxiliarySource.getValue() ||
	(!hasCompactPayload &&
	 (this->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	  this->needsRealization())) ||
	source_has_auxiliary_children(this)) {
	this->d->compiledCompactStructureSignature = 0;
	this->d->compiledCompactStyleSignature = 0;
	this->d->compiledCompactSemanticSignature = 0;
	this->d->compiledCompactHiddenSignature = 0;
	this->d->compiledCompactSelectedSignature = 0;
	this->d->compiledCompactUnpickableSignature = 0;
	this->d->compiledAssemblyNodeId = sourceNodeId;
	this->d->compiledAssemblyDirty = FALSE;
	return 0;
    }

    if (!this->d->compiledAssembly) {
	this->d->compiledAssembly = new SoBRLCadAssembly;
	this->d->compiledAssembly->ref();
    }

    const uint64_t structureSignature = hasCompactPayload ?
	compact_structure_signature(this->d->compactIndex) : 0;
    const uint64_t semanticSignature = hasCompactPayload ?
	compact_semantic_signature(this->d->compactIndex) : 0;
    const uint64_t styleSignature = hasCompactPayload ?
	compact_style_signature(this->d->compactIndex) : 0;
    if (hasCompactPayload && this->d->compiledCompactStructureSignature != 0 &&
	this->d->compiledCompactStructureSignature == structureSignature &&
	this->d->compiledAssembly->instanceCount() ==
	    this->d->compactIndex->instances.size() &&
	this->d->compiledAssembly->partCount() == this->d->compactIndex->parts.size()) {
	if (this->d->compiledCompactStyleSignature != styleSignature) {
	    std::vector<Obol::InstanceStyleUpdate> styles;
	    styles.reserve(this->d->compactIndex->instances.size());
	    for (const Obol::InstanceUpdate &instance :
		 this->d->compactIndex->instances) {
		Obol::InstanceStyleUpdate style;
		style.instance = instance.instance;
		style.style = instance.record.style;
		styles.push_back(style);
	    }
	    this->d->compiledAssembly->updateInstanceStyles(styles);
	    this->d->compiledCompactStyleSignature = styleSignature;
	}
	const uint64_t hiddenSignature = compact_state_signature(
	    this->d->compactIndex->hiddenInstances);
	if (this->d->compiledCompactHiddenSignature != hiddenSignature) {
	    this->d->compiledAssembly->setHiddenInstances(
		this->d->compactIndex->hiddenInstances);
	    this->d->compiledCompactHiddenSignature = hiddenSignature;
	}
	const uint64_t selectedSignature = compact_state_signature(
	    this->d->compactIndex->selectedInstances);
	if (this->d->compiledCompactSelectedSignature != selectedSignature) {
	    this->d->compiledAssembly->setSelectedInstances(
		this->d->compactIndex->selectedInstances);
	    this->d->compiledCompactSelectedSignature = selectedSignature;
	}
	const uint64_t unpickableSignature = compact_state_signature(
	    this->d->compactIndex->unpickableInstances);
	if (this->d->compiledCompactUnpickableSignature != unpickableSignature) {
	    this->d->compiledAssembly->setUnpickableInstances(
		this->d->compactIndex->unpickableInstances);
	    this->d->compiledCompactUnpickableSignature = unpickableSignature;
	}
	if (this->d->compiledCompactSemanticSignature != semanticSignature) {
	    for (const BObolCompactInstanceEntry &entry :
		 this->d->compactIndex->entries) {
		SoBRLCadAssembly::InstanceSemantic semantic = entry.semantic;
		semantic.sourceInstanceKey = compact_instance_identity(entry);
		this->d->compiledAssembly->setInstanceSemantic(entry.instance, semantic);
	    }
	    this->d->compiledCompactSemanticSignature = semanticSignature;
	}
	compact_assembly_draw_mode(this->d->compiledAssembly, this,
	    this->d->compactIndex);
	this->d->compiledAssemblyActive = TRUE;
	this->ensureCompiledAssemblyChild();
	this->d->compiledAssemblyNodeId = this->getNodeId();
	this->d->compiledAssemblyDirty = FALSE;
	return 1;
    }

    this->d->compiledAssembly->beginUpdate();
    this->d->compiledAssembly->clear();
    this->d->compiledAssembly->clearSemanticMap();

    if (this->d->compactIndexActive && this->d->compactIndex &&
	!this->d->compactIndex->instances.empty()) {
	std::vector<Obol::SharedPartUpdate> parts;
	parts.reserve(this->d->compactIndex->parts.size());
	for (const BObolCompactPartReference &partRef :
	     this->d->compactIndex->parts) {
	    if (!partRef.geometry)
		continue;
	    Obol::SharedPartUpdate part;
	    part.part = partRef.part;
	    part.geometry = partRef.geometry;
	    parts.push_back(part);
	}
	this->d->compiledAssembly->upsertSharedParts(parts);
	this->d->compiledAssembly->upsertInstances(this->d->compactIndex->instances);
	this->d->compiledAssembly->setHiddenInstances(
	    this->d->compactIndex->hiddenInstances);
	this->d->compiledAssembly->setSelectedInstances(
	    this->d->compactIndex->selectedInstances);
	this->d->compiledAssembly->setUnpickableInstances(
	    this->d->compactIndex->unpickableInstances);
	this->d->compiledCompactHiddenSignature = compact_state_signature(
	    this->d->compactIndex->hiddenInstances);
	this->d->compiledCompactSelectedSignature = compact_state_signature(
	    this->d->compactIndex->selectedInstances);
	this->d->compiledCompactUnpickableSignature = compact_state_signature(
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
	this->d->compiledCompactStructureSignature = structureSignature;
	this->d->compiledCompactStyleSignature = styleSignature;
	this->d->compiledCompactSemanticSignature = semanticSignature;
	this->d->compiledAssembly->endUpdate();
	this->ensureCompiledAssemblyChild();
	this->d->compiledAssemblyNodeId = this->getNodeId();
	this->d->compiledAssemblyDirty = FALSE;
	return 1;
    }

    /* The retained compact signature has no meaning for a child-graph build. */
    this->d->compiledCompactStructureSignature = 0;
    this->d->compiledCompactStyleSignature = 0;
    this->d->compiledCompactSemanticSignature = 0;
    this->d->compiledCompactHiddenSignature = 0;
    this->d->compiledCompactSelectedSignature = 0;
    this->d->compiledCompactUnpickableSignature = 0;
    cad_build_data data;
    data.source = this;
    data.assembly = this->d->compiledAssembly;
    data.ordinal = 0;
    data.unsupported = 0;
    data.wireCount = 0;
    data.shadedCount = 0;
    data.storeSemantics = 1;

    const SbMatrix identity = SbMatrix::identity();
    for (int i = 0; i < this->getNumChildren() && !data.unsupported; i++)
	cad_collect_realized_node(data, this->getChild(i), identity);

    if (!data.unsupported && !data.instances.empty()) {
	this->d->compiledAssembly->upsertParts(data.parts);
	this->d->compiledAssembly->upsertInstances(data.instances);
	this->d->compiledAssembly->setHiddenInstances(data.hiddenInstances);
	this->d->compiledAssembly->setSelectedInstances(data.selectedInstances);
	this->d->compiledAssembly->setUnpickableInstances(
	    data.unpickableInstances);
	if (source_record_draw_mode(this) == BOBOL_LOD_DRAW_HIDDEN_LINE)
	    this->d->compiledAssembly->drawMode = SoCADAssembly::HIDDEN_LINE;
	else if (data.shadedCount > 0 && data.wireCount > 0)
	    this->d->compiledAssembly->drawMode = SoCADAssembly::SHADED_WITH_EDGES;
	else if (data.shadedCount > 0)
	    this->d->compiledAssembly->drawMode = SoCADAssembly::SHADED;
	else
	    this->d->compiledAssembly->drawMode = SoCADAssembly::WIREFRAME;
	this->d->compiledAssemblyActive = TRUE;
    }

    this->d->compiledAssembly->endUpdate();
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
    if (!state || !state->assembly || !this->visible.getValue() ||
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
		Obol::SharedPartUpdate part;
		part.part = partRef.part;
		part.geometry = partRef.geometry;
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
		state->assembly->setInstanceSemantic(entry.instance, semantic);
	    }
	}
	state->wireCount += this->d->compactIndex->wireCount;
	state->shadedCount += this->d->compactIndex->shadedCount;
	return 1;
    }

    cad_build_data data;
    data.source = this;
    data.assembly = state->assembly;
    data.ordinal = 0;
    data.unsupported = 0;
    data.wireCount = 0;
    data.shadedCount = 0;
    data.storeSemantics = 0;

    const SbMatrix identity = SbMatrix::identity();
    for (int i = 0; i < this->getNumChildren() && !data.unsupported; i++)
	cad_collect_realized_node(data, this->getChild(i), identity);
    if (data.unsupported || data.instances.empty())
	return 0;

    for (Obol::PartUpdate &part : data.parts) {
	if (state->partIds.insert(part.part).second) {
	    Obol::SharedPartUpdate sharedPart;
	    sharedPart.part = part.part;
	    sharedPart.geometry =
		std::make_shared<const Obol::PartGeometry>(std::move(part.geometry));
	    state->parts.push_back(std::move(sharedPart));
	}
    }
    state->instances.insert(state->instances.end(), data.instances.begin(),
	data.instances.end());
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

    this->attachFieldSensors();
}

SoBRLDatabaseSource::~SoBRLDatabaseSource(void)
{
    this->detachFieldSensors();
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

void
SoBRLDatabaseSource::fieldSensorCB(void *data, SoSensor *sensor)
{
    SoBRLDatabaseSource *source = static_cast<SoBRLDatabaseSource *>(data);
    if (!source)
	return;

    if (sensor == source->d->inputsRevisionSensor)
	source->markStale(STALE_INPUTS);
    else if (sensor == source->d->viewRevisionSensor ||
	     sensor == source->d->lodBotThresholdSensor)
	source->markStale(STALE_VIEW);
    else if (sensor == source->d->drawModeSensor ||
	     sensor == source->d->representationModeSensor)
	source->markStale(STALE_DRAW);
    else if (sensor == source->d->tessellationAbsTolSensor ||
	     sensor == source->d->tessellationRelTolSensor ||
	     sensor == source->d->tessellationNormTolSensor)
	source->markStale(STALE_TESSELLATION);
    else
	source->markStale(STALE_SOURCE);
}

void
SoBRLDatabaseSource::attachFieldSensors(void)
{
    this->d->pathSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->d->pathSensor->setPriority(0);
    this->d->pathSensor->attach(&this->path);

    this->d->instanceKeySensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->d->instanceKeySensor->setPriority(0);
    this->d->instanceKeySensor->attach(&this->instanceKey);

    this->d->representationKeySensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->d->representationKeySensor->setPriority(0);
    this->d->representationKeySensor->attach(&this->representationKey);

    this->d->representationModeSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->d->representationModeSensor->setPriority(0);
    this->d->representationModeSensor->attach(&this->representationMode);

    this->d->drawModeSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->d->drawModeSensor->setPriority(0);
    this->d->drawModeSensor->attach(&this->drawMode);

    this->d->tessellationAbsTolSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->d->tessellationAbsTolSensor->setPriority(0);
    this->d->tessellationAbsTolSensor->attach(&this->tessellationAbsTol);

    this->d->tessellationRelTolSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->d->tessellationRelTolSensor->setPriority(0);
    this->d->tessellationRelTolSensor->attach(&this->tessellationRelTol);

    this->d->tessellationNormTolSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->d->tessellationNormTolSensor->setPriority(0);
    this->d->tessellationNormTolSensor->attach(&this->tessellationNormTol);

    this->d->lodBotThresholdSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->d->lodBotThresholdSensor->setPriority(0);
    this->d->lodBotThresholdSensor->attach(&this->lodBotThreshold);

    this->d->sourceRevisionSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->d->sourceRevisionSensor->setPriority(0);
    this->d->sourceRevisionSensor->attach(&this->sourceRevision);

    this->d->inputsRevisionSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->d->inputsRevisionSensor->setPriority(0);
    this->d->inputsRevisionSensor->attach(&this->inputsRevision);

    this->d->viewRevisionSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->d->viewRevisionSensor->setPriority(0);
    this->d->viewRevisionSensor->attach(&this->viewRevision);
}

void
SoBRLDatabaseSource::detachFieldSensors(void)
{
    if (this->d->pathSensor)
	this->d->pathSensor->detach();
    delete this->d->pathSensor;
    this->d->pathSensor = NULL;
    if (this->d->instanceKeySensor)
	this->d->instanceKeySensor->detach();
    delete this->d->instanceKeySensor;
    this->d->instanceKeySensor = NULL;
    if (this->d->representationKeySensor)
	this->d->representationKeySensor->detach();
    delete this->d->representationKeySensor;
    this->d->representationKeySensor = NULL;
    if (this->d->representationModeSensor)
	this->d->representationModeSensor->detach();
    delete this->d->representationModeSensor;
    this->d->representationModeSensor = NULL;
    if (this->d->drawModeSensor)
	this->d->drawModeSensor->detach();
    delete this->d->drawModeSensor;
    this->d->drawModeSensor = NULL;
    if (this->d->tessellationAbsTolSensor)
	this->d->tessellationAbsTolSensor->detach();
    delete this->d->tessellationAbsTolSensor;
    this->d->tessellationAbsTolSensor = NULL;
    if (this->d->tessellationRelTolSensor)
	this->d->tessellationRelTolSensor->detach();
    delete this->d->tessellationRelTolSensor;
    this->d->tessellationRelTolSensor = NULL;
    if (this->d->tessellationNormTolSensor)
	this->d->tessellationNormTolSensor->detach();
    delete this->d->tessellationNormTolSensor;
    this->d->tessellationNormTolSensor = NULL;
    if (this->d->lodBotThresholdSensor)
	this->d->lodBotThresholdSensor->detach();
    delete this->d->lodBotThresholdSensor;
    this->d->lodBotThresholdSensor = NULL;
    if (this->d->sourceRevisionSensor)
	this->d->sourceRevisionSensor->detach();
    delete this->d->sourceRevisionSensor;
    this->d->sourceRevisionSensor = NULL;
    if (this->d->inputsRevisionSensor)
	this->d->inputsRevisionSensor->detach();
    delete this->d->inputsRevisionSensor;
    this->d->inputsRevisionSensor = NULL;
    if (this->d->viewRevisionSensor)
	this->d->viewRevisionSensor->detach();
    delete this->d->viewRevisionSensor;
    this->d->viewRevisionSensor = NULL;
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
    this->d->compiledCompactStructureSignature = 0;
    this->d->compiledCompactStyleSignature = 0;
    this->d->compiledCompactSemanticSignature = 0;
    this->d->compiledCompactHiddenSignature = 0;
    this->d->compiledCompactSelectedSignature = 0;
    this->d->compiledCompactUnpickableSignature = 0;
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
    if (++this->d->cadBatchRevision == 0)
	this->d->cadBatchRevision = 1;
    this->d->cadBatchDeltas.clear();
    this->d->cadBatchDeltaEntryCount = 0;
    this->d->cadBatchDeltaFloorRevision =
	this->d->cadBatchRevision;
}

void
SoBRLDatabaseSource::markCadBatchDirty(
    const std::vector<size_t> &entryIndices)
{
    if (entryIndices.empty())
	return;
    if (++this->d->cadBatchRevision == 0)
	this->d->cadBatchRevision = 1;

    Impl::CadBatchDelta delta;
    delta.revision = this->d->cadBatchRevision;
    delta.entryIndices = entryIndices;
    std::sort(delta.entryIndices.begin(), delta.entryIndices.end());
    delta.entryIndices.erase(std::unique(delta.entryIndices.begin(),
	delta.entryIndices.end()), delta.entryIndices.end());
    this->d->cadBatchDeltaEntryCount += delta.entryIndices.size();
    this->d->cadBatchDeltas.push_back(std::move(delta));

    /* Multiple views may acknowledge the same source independently.  Keep
     * the journal non-consuming and bounded; a view that falls behind the
     * retained floor performs one authoritative scan. */
    static const size_t maxDeltaBatches = 256;
    static const size_t maxDeltaEntries = 65536;
    while (this->d->cadBatchDeltas.size() > maxDeltaBatches ||
	this->d->cadBatchDeltaEntryCount > maxDeltaEntries) {
	const Impl::CadBatchDelta &front =
	    this->d->cadBatchDeltas.front();
	this->d->cadBatchDeltaEntryCount -= front.entryIndices.size();
	this->d->cadBatchDeltaFloorRevision =
	    std::max(this->d->cadBatchDeltaFloorRevision,
		front.revision);
	this->d->cadBatchDeltas.pop_front();
    }
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
    if (++this->d->displayMeshLodRevision == 0)
	this->d->displayMeshLodRevision = 1;
    this->d->displayMeshLodDeltas.clear();
    this->d->displayMeshLodDeltaEntryCount = 0;
    this->d->displayMeshLodDeltaFloorRevision =
	this->d->displayMeshLodRevision;
}

void
SoBRLDatabaseSource::markDisplayMeshLodDirty(
    const std::vector<size_t> &entryIndices)
{
    if (entryIndices.empty())
	return;
    if (++this->d->displayMeshLodRevision == 0)
	this->d->displayMeshLodRevision = 1;

    Impl::DisplayMeshLodDelta delta;
    delta.revision = this->d->displayMeshLodRevision;
    delta.entryIndices = entryIndices;
    std::sort(delta.entryIndices.begin(), delta.entryIndices.end());
    delta.entryIndices.erase(std::unique(delta.entryIndices.begin(),
	delta.entryIndices.end()), delta.entryIndices.end());
    this->d->displayMeshLodDeltaEntryCount += delta.entryIndices.size();
    this->d->displayMeshLodDeltas.push_back(std::move(delta));

    /* A visibility/edit burst normally records a handful of entries.  Bound
     * pathological history without penalizing independent views: a lagging
     * view receives FALSE below and performs one authoritative rescan. */
    static const size_t maxDeltaBatches = 256;
    static const size_t maxDeltaEntries = 65536;
    while (this->d->displayMeshLodDeltas.size() > maxDeltaBatches ||
	this->d->displayMeshLodDeltaEntryCount > maxDeltaEntries) {
	const Impl::DisplayMeshLodDelta &front =
	    this->d->displayMeshLodDeltas.front();
	this->d->displayMeshLodDeltaEntryCount -=
	    front.entryIndices.size();
	this->d->displayMeshLodDeltaFloorRevision =
	    std::max(this->d->displayMeshLodDeltaFloorRevision,
		front.revision);
	this->d->displayMeshLodDeltas.pop_front();
    }
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
    this->d->compactExpectedInstanceCount = 0;
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
    this->d->compactExpectedInstanceCount = 0;
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
    BObolCompactInstanceIndex *index, SbBool occurrenceRegistry)
{
    if (!index)
	return;
    const BObolCompactInstanceIndex *previous = this->d->compactIndex ?
	this->d->compactIndex : this->d->previousCompactIndex;
    compact_merge_runtime_state(previous, index);
    this->clearCompactInstanceIndex();
    this->d->compactIndex = index;
    this->d->compactIndexActive = TRUE;
    this->d->compactOccurrenceRegistry = occurrenceRegistry;
    this->d->compactExpectedInstanceCount = index->entries.size();
    this->d->displayMeshLodContractRevisionValid =
	index->sourceMeshRequestCount > 0 ? TRUE : FALSE;
    this->d->displayMeshLodContractSourceRevision =
	this->sourceRevision.getValue();
    this->d->displayMeshLodContractInputsRevision =
	this->inputsRevision.getValue();
    if (getenv("BOBOL_LOD_TRACE_SOURCE_CONTRACT"))
	bu_log("BObol LoD source contract installed compact index path=%s "
	       "entries=%zu requests=%zu registry=%d draw=%d threshold=%u "
	       "view_dependent=%d mesh_lod=%d\n",
	       this->path.getValue().getString(), index->entries.size(),
	       index->sourceMeshRequestCount,
	       occurrenceRegistry ? 1 : 0,
	       source_record_draw_mode(this),
	       this->lodBotThreshold.getValue(),
	       this->realizationViewDependent.getValue() ? 1 : 0,
	       this->realizationMeshLodEnabled.getValue() ? 1 : 0);
    (void)this->reapplyCompactInstanceVisibilityFrontier();
    (void)this->reapplyCompactInstanceSelectedPaths();

    delete this->d->previousCompactIndex;
    this->d->previousCompactIndex = NULL;
}

void
SoBRLDatabaseSource::markStale(void)
{
    this->markStale(STALE_SOURCE);
}

void
SoBRLDatabaseSource::markStale(uint32_t reason)
{
    if (!reason)
	return;
    const SbBool hadCurrentDisplayMeshLodContract =
	this->hasDisplayMeshLodRequests();
    this->stale = TRUE;
    this->staleReason = this->staleReason.getValue() | reason;
    /*
     * A compact source-mesh request identifies immutable database geometry;
     * the camera, draw channel, tessellation display policy, and PoP cut are
     * deliberately supplied later by the view-local LoD planner.  Invalidating
     * that contract on STALE_VIEW made a camera update race the final streamed
     * adoption: the source became current again, but wireframe sources (which
     * have no shaded fallback mesh) disappeared from the planner and left its
     * submission cursor spinning forever.
     *
     * Only an identity/input/database change can make the source request
     * itself unsafe.  Other stale reasons may replace presentation data in
     * the background while the existing request and retained PoP generations
     * remain valid and useful.
     */
    const uint32_t sourceContractReasons =
	STALE_SOURCE | STALE_INPUTS | STALE_DATABASE;
    if (reason & sourceContractReasons)
	this->d->displayMeshLodContractRevisionValid = FALSE;
    if (getenv("BOBOL_LOD_TRACE_SOURCE_CONTRACT") &&
	hadCurrentDisplayMeshLodContract &&
	!this->hasDisplayMeshLodRequests())
	bu_log("BObol LoD source contract invalidated path=%s reason=%u "
	       "source_revision=%u inputs_revision=%u requests=%zu\n",
	       this->path.getValue().getString(), reason,
	       this->sourceRevision.getValue(),
	       this->inputsRevision.getValue(),
	       this->d->compactIndex ?
		   this->d->compactIndex->sourceMeshRequestCount : 0);
    this->realizationStatus = UNREALIZED;
    this->realizationDiagnostic = "";
    /* Retain immutable compact geometry until its replacement is ready.  A
     * stale source describes realization currency, not display validity. */
    this->markCompiledAssemblyDirty();
    this->markCadBatchDirty();
    this->syncRealizedShapeOwnerState();
}

int
SoBRLDatabaseSource::setDrawModeState(int nextDrawMode)
{
    if (nextDrawMode != SHADED)
	nextDrawMode = WIREFRAME;
    if (this->drawMode.getValue() == nextDrawMode)
	return 0;

    const SbBool preserveExternalRealization =
	(this->realizationRoleFlags.getValue() & REALIZATION_ROLE_EXTERNAL) &&
	this->realizationStatus.getValue() == REALIZED &&
	!this->stale.getValue() &&
	((nextDrawMode == SHADED && this->hasRealizedMeshGeometry()) ||
	 (nextDrawMode == WIREFRAME && this->hasRealizedWireGeometry()));

    if (this->d->drawModeSensor)
	this->d->drawModeSensor->detach();
    this->drawMode = nextDrawMode;
    if (this->d->drawModeSensor)
	this->d->drawModeSensor->attach(&this->drawMode);
    if (!preserveExternalRealization)
	this->markStale(STALE_DRAW);
    this->syncRealizedShapeOwnerState();
    return 1;
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

int
SoBRLDatabaseSource::setDisplayNameState(const char *name)
{
    const char *nextName = name ? name : "";
    if (database_source_string_equal(this->displayName.getValue(), nextName))
	return 0;

    this->displayName = nextName;
    this->syncRealizedShapeOwnerState();
    return 1;
}

int
SoBRLDatabaseSource::setRepresentationState(
    const char *sourceRepresentationKey,
    int sourceRepresentationMode)
{
    const char *nextKey = sourceRepresentationKey ?
			  sourceRepresentationKey : "";
    const int nextMode = sourceRepresentationMode;
    if (database_source_string_equal(this->representationKey.getValue(),
				     nextKey) && this->representationMode.getValue() == nextMode)
	return 0;

    if (this->d->representationKeySensor)
	this->d->representationKeySensor->detach();
    if (this->d->representationModeSensor)
	this->d->representationModeSensor->detach();
    this->representationKey = nextKey;
    this->representationMode = nextMode;
    if (this->d->representationKeySensor)
	this->d->representationKeySensor->attach(&this->representationKey);
    if (this->d->representationModeSensor)
	this->d->representationModeSensor->attach(&this->representationMode);

    this->markStale(STALE_DRAW);
    const char *sourcePath = this->path.getValue().getString();
    for (int i = 0; i < this->getNumChildren(); i++)
	retarget_realized_node_source(this->getChild(i), this,
				      sourcePath, sourcePath,
				      this->sourceRevision.getValue());
    this->realizationIdentity = source_realization_identity(this);
    return 1;
}

int
SoBRLDatabaseSource::setHierarchyState(
    const char *sourceParentInstanceKey,
    uint32_t sourceOccurrenceIndex,
    int sourceBooleanOperation)
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

    this->parentInstanceKey = nextParent;
    this->occurrenceIndex = sourceOccurrenceIndex;
    this->booleanOperation = nextOperation;
    this->markCompiledAssemblyDirty();
    this->markCadBatchDirty();
    return 1;
}

int
SoBRLDatabaseSource::setMaterialPolicyState(int nextMaterialPolicy)
{
    if (nextMaterialPolicy != MATERIAL_DATABASE)
	nextMaterialPolicy = MATERIAL_INHERIT;
    if (this->materialPolicy.getValue() == nextMaterialPolicy)
	return 0;

    this->materialPolicy = nextMaterialPolicy;
    this->markCadBatchDirty();
    this->syncRealizedShapeOwnerState();
    return 1;
}

int
SoBRLDatabaseSource::setRealizationState(int nextStatus,
	uint32_t nextRealizedSourceRevision,
	uint32_t nextRealizedInputsRevision,
	uint32_t nextStaleReason,
	const char *diagnostic)
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
    SbString nextDiagnostic = diagnostic ? diagnostic : "";

    int changed = 0;
    if (this->realizationStatus.getValue() != nextStatus) {
	this->realizationStatus = nextStatus;
	changed = 1;
    }
    if (this->stale.getValue() != nextStale) {
	this->stale = nextStale;
	changed = 1;
    }
    if (this->staleReason.getValue() != nextStaleReason) {
	this->staleReason = nextStaleReason;
	changed = 1;
    }
    if (bu_strcmp(this->realizationDiagnostic.getValue().getString(),
	       nextDiagnostic.getString()) != 0) {
	this->realizationDiagnostic = nextDiagnostic;
	changed = 1;
    }
    if (realized &&
	this->realizedSourceRevision.getValue() !=
	nextRealizedSourceRevision) {
	this->realizedSourceRevision = nextRealizedSourceRevision;
	changed = 1;
    }
    if (realized &&
	this->realizedInputsRevision.getValue() !=
	nextRealizedInputsRevision) {
	this->realizedInputsRevision = nextRealizedInputsRevision;
	changed = 1;
    }
    if (realized &&
	this->realizedViewRevision.getValue() !=
	this->viewRevision.getValue()) {
	this->realizedViewRevision = this->viewRevision.getValue();
	changed = 1;
    }
    if (realized &&
	this->realizedRevision.getValue() != nextRealizedSourceRevision) {
	this->realizedRevision = nextRealizedSourceRevision;
	changed = 1;
    }
    const SbBool isBatchEligible = realized && !nextStale;
    if (wasBatchEligible != isBatchEligible)
	this->markCadBatchDirty();
    if (changed)
	this->syncRealizedShapeOwnerState();
    return changed;
}

int
SoBRLDatabaseSource::setRealizationRoleFlags(int roleFlags)
{
    const int validFlags = REALIZATION_ROLE_CSG | REALIZATION_ROLE_MESH |
			   REALIZATION_ROLE_EXTERNAL;
    roleFlags &= validFlags;
    if (this->realizationRoleFlags.getValue() == roleFlags)
	return 0;

    this->realizationRoleFlags = roleFlags;
    this->syncRealizedShapeOwnerState();
    return 1;
}

int
SoBRLDatabaseSource::setRealizationViewPolicy(SbBool viewDependent,
	SbBool csgLodEnabled,
	SbBool meshLodEnabled,
	float viewScale,
	float lodScale,
	int viewWidth,
	int viewHeight,
	uint32_t botThreshold,
	float curveScale,
	float pointScale)
{
    int changed = 0;
    int activePolicyChanged = 0;
    if (this->realizationViewDependent.getValue() != viewDependent) {
	this->realizationViewDependent = viewDependent;
	activePolicyChanged = 1;
	changed = 1;
    }
    if (this->realizationCsgLodEnabled.getValue() != csgLodEnabled) {
	this->realizationCsgLodEnabled = csgLodEnabled;
	activePolicyChanged = 1;
	changed = 1;
    }
    if (this->realizationMeshLodEnabled.getValue() != meshLodEnabled) {
	this->realizationMeshLodEnabled = meshLodEnabled;
	activePolicyChanged = 1;
	changed = 1;
    }
    if (database_source_float_different(
	    this->realizationViewScale.getValue(), viewScale)) {
	this->realizationViewScale = viewScale;
	activePolicyChanged = 1;
	changed = 1;
    }
    if (lodScale <= 0.0f)
	lodScale = 1.0f;
    if (database_source_float_different(
	    this->realizationLodScale.getValue(), lodScale)) {
	this->realizationLodScale = lodScale;
	activePolicyChanged = 1;
	changed = 1;
    }
    if (this->realizationViewWidth.getValue() != viewWidth) {
	this->realizationViewWidth = viewWidth;
	activePolicyChanged = 1;
	changed = 1;
    }
    if (this->realizationViewHeight.getValue() != viewHeight) {
	this->realizationViewHeight = viewHeight;
	activePolicyChanged = 1;
	changed = 1;
    }
    if (this->realizationBotThreshold.getValue() != botThreshold) {
	this->realizationBotThreshold = botThreshold;
	activePolicyChanged = 1;
	changed = 1;
    }
    if (this->lodBotThreshold.getValue() != botThreshold) {
	if (this->d->lodBotThresholdSensor)
	    this->d->lodBotThresholdSensor->detach();
	this->lodBotThreshold = botThreshold;
	if (this->d->lodBotThresholdSensor)
	    this->d->lodBotThresholdSensor->attach(&this->lodBotThreshold);
	activePolicyChanged = 1;
	changed = 1;
    }
    if (database_source_float_different(
	    this->realizationCurveScale.getValue(), curveScale)) {
	this->realizationCurveScale = curveScale;
	activePolicyChanged = 1;
	changed = 1;
    }
    if (database_source_float_different(
	    this->realizationPointScale.getValue(), pointScale)) {
	this->realizationPointScale = pointScale;
	activePolicyChanged = 1;
	changed = 1;
    }
    if (activePolicyChanged)
	this->markStale(STALE_VIEW);
    if (changed)
	this->syncRealizedShapeOwnerState();
    return changed;
}

int
SoBRLDatabaseSource::setDatabaseMetadataState(SbBool metadataValid,
	int regionId,
	int airCode,
	int materialId,
	int los,
	SbBool metadataMaterialColorValid,
	const SbColor &metadataMaterialColor,
	const SbString &metadataMaterialShader)
{
    int changed = 0;
    if (this->databaseMetadataValid.getValue() != metadataValid) {
	this->databaseMetadataValid = metadataValid;
	changed = 1;
    }
    if (this->databaseRegionId.getValue() != regionId) {
	this->databaseRegionId = regionId;
	changed = 1;
    }
    if (this->databaseAirCode.getValue() != airCode) {
	this->databaseAirCode = airCode;
	changed = 1;
    }
    if (this->databaseMaterialId.getValue() != materialId) {
	this->databaseMaterialId = materialId;
	changed = 1;
    }
    if (this->databaseLos.getValue() != los) {
	this->databaseLos = los;
	changed = 1;
    }
    if (this->databaseMaterialColorValid.getValue() !=
	metadataMaterialColorValid) {
	this->databaseMaterialColorValid = metadataMaterialColorValid;
	changed = 1;
    }
    if (metadataMaterialColorValid &&
	!database_source_color_equal(
	    this->databaseMaterialColor.getValue(), metadataMaterialColor)) {
	this->databaseMaterialColor = metadataMaterialColor;
	changed = 1;
    }
    if (!metadataMaterialColorValid &&
	!database_source_color_equal(this->databaseMaterialColor.getValue(),
				     SbColor(1.0f, 1.0f, 1.0f))) {
	this->databaseMaterialColor = SbColor(1.0f, 1.0f, 1.0f);
	changed = 1;
    }
    if (bu_strcmp(this->databaseMaterialShader.getValue().getString(),
	       metadataMaterialShader.getString()) != 0) {
	this->databaseMaterialShader = metadataMaterialShader;
	changed = 1;
    }

    if (changed) {
	this->markCadBatchDirty();
	this->syncRealizedShapeOwnerState();
    }
    return changed;
}

int
SoBRLDatabaseSource::refreshMaterialColorFromDatabase(
    uint32_t nextMaterialRevision,
    struct db_i *overrideDbip)
{
    /* A material revision is the synchronization stamp for the effective
     * display color.  The primary scene is refreshed before attached-scene
     * observers run, so avoid repeating the comparatively expensive full-path
     * database lookup when an observer targets that same scene. */
    if (this->materialRevision.getValue() == nextMaterialRevision)
	return 0;

    struct db_i *colorDbip = overrideDbip ? overrideDbip : this->d->dbip;
    if (this->d->compactIndex) {
	if (!colorDbip)
	    return 0;
	BObolMaterialColorSweep sweep(colorDbip);
	int changed = 0;
	BObolMaterialPathState sourceResolved;
	if (sweep.resolve(this->path.getValue().getString(), sourceResolved) &&
	    (!this->materialColorValid.getValue() ||
	     !database_source_color_equal(this->materialColor.getValue(),
		 sourceResolved.color))) {
	    /* The aggregate summary tracks its own effective color, but compact
	     * entries retain their individual full-path material colors below. */
	    this->materialColorValid = TRUE;
	    this->materialColor = sourceResolved.color;
	    changed = 1;
	}
	for (BObolCompactInstanceEntry &entry : this->d->compactIndex->entries) {
	    BObolMaterialPathState resolved;
	    if (!sweep.resolve(entry.semantic.path.getString(), resolved))
		continue;
	    const SbString shader(resolved.shader.c_str());
	    if (entry.semantic.regionId == resolved.regionId &&
		entry.semantic.airCode == resolved.airCode &&
		entry.semantic.materialId == resolved.materialId &&
		entry.semantic.los == resolved.los &&
		entry.semantic.materialColorValid &&
		database_source_color_equal(entry.semantic.materialColor,
		    resolved.color) &&
		bu_strcmp(entry.semantic.materialShader.getString(),
		    shader.getString()) == 0)
		continue;
	    entry.semantic.regionId = resolved.regionId;
	    entry.semantic.airCode = resolved.airCode;
	    entry.semantic.materialId = resolved.materialId;
	    entry.semantic.los = resolved.los;
	    entry.semantic.materialColorValid = TRUE;
	    entry.semantic.materialColor = resolved.color;
	    entry.semantic.materialShader = shader;
	    entry.appearanceRevision = compact_next_revision(
		entry.appearanceRevision);
	    changed = 1;
	}
	if (!changed) {
	    this->materialRevision = nextMaterialRevision;
	    return 0;
	}
	this->materialRevision = nextMaterialRevision;
	/* The sweep above updates occurrence semantics.  Rebuild each retained
	 * instance's effective style from those new colors before publishing the
	 * compact batch; otherwise the summaries report the new table while the
	 * renderer continues using the previous compiled color. */
	for (BObolCompactInstanceEntry &entry : this->d->compactIndex->entries)
	    compact_sync_entry_from_source(entry, this);
	this->rebuildCompactInstanceDisplayState(FALSE);
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
	this->touch();
	return changed;
    }
    SbColor nextMaterialColor(1.0f, 1.0f, 1.0f);
    if (!bobol_database_source_path_material_color(
	    colorDbip,
	    this->path.getValue().getString(),
	    nextMaterialColor))
	return 0;

    return this->setDisplayState(
	FALSE,
	this->sourceRevision.getValue(),
	this->inputsRevision.getValue(),
	this->visible.getValue(),
	this->selected.getValue(),
	this->highlighted.getValue(),
	this->lineStyle.getValue(),
	this->lineWidth.getValue(),
	this->transparency.getValue(),
	this->colorOverride.getValue(),
	this->color.getValue(),
	TRUE,
	nextMaterialColor,
	nextMaterialRevision);
}

int
SoBRLDatabaseSource::setDisplayState(SbBool sourceRevisionValid,
				     uint32_t nextSourceRevision,
				     uint32_t nextInputsRevision,
				     SbBool nextVisible,
				     SbBool nextSelected,
				     SbBool nextHighlighted,
				     int nextLineStyle,
				     int nextLineWidth,
				     float nextTransparency,
				     SbBool nextColorOverride,
				     const SbColor &nextColor,
				     SbBool nextMaterialColorValid,
				     const SbColor &nextMaterialColor,
				     uint32_t nextMaterialRevision)
{
    int changed = 0;
    int appearanceChanged = 0;
    int visibilityChanged = 0;
    int selectedChanged = 0;
    int highlightedChanged = 0;
    const SbBool hadMaterialOverride = this->materialColorValid.getValue();
    if (sourceRevisionValid &&
	this->sourceRevision.getValue() != nextSourceRevision) {
	this->sourceRevision = nextSourceRevision;
	changed = 1;
    }
    if (this->inputsRevision.getValue() != nextInputsRevision) {
	this->inputsRevision = nextInputsRevision;
	changed = 1;
    }
    if (this->visible.getValue() != nextVisible) {
	this->visible = nextVisible;
	visibilityChanged = 1;
	changed = 1;
    }
    if (this->selected.getValue() != nextSelected) {
	this->selected = nextSelected;
	selectedChanged = 1;
	changed = 1;
    }
    if (this->highlighted.getValue() != nextHighlighted) {
	this->highlighted = nextHighlighted;
	highlightedChanged = 1;
	changed = 1;
    }
    if (this->lineStyle.getValue() != nextLineStyle) {
	this->lineStyle = nextLineStyle;
	appearanceChanged = 1;
	changed = 1;
    }
    if (this->lineWidth.getValue() != nextLineWidth) {
	this->lineWidth = nextLineWidth;
	appearanceChanged = 1;
	changed = 1;
    }
    if (database_source_float_different(this->transparency.getValue(),
					nextTransparency)) {
	this->transparency = nextTransparency;
	appearanceChanged = 1;
	changed = 1;
    }
    if (this->colorOverride.getValue() != nextColorOverride) {
	this->colorOverride = nextColorOverride;
	appearanceChanged = 1;
	changed = 1;
    }
    if (nextColorOverride &&
	!database_source_color_equal(this->color.getValue(), nextColor)) {
	this->color = nextColor;
	appearanceChanged = 1;
	changed = 1;
    }
    if (this->materialColorValid.getValue() != nextMaterialColorValid) {
	this->materialColorValid = nextMaterialColorValid;
	appearanceChanged = 1;
	changed = 1;
    }
    if (nextMaterialColorValid &&
	!database_source_color_equal(this->materialColor.getValue(),
				     nextMaterialColor)) {
	this->materialColor = nextMaterialColor;
	appearanceChanged = 1;
	changed = 1;
    }
    if (this->materialRevision.getValue() != nextMaterialRevision) {
	this->materialRevision = nextMaterialRevision;
	changed = 1;
    }
    if (hadMaterialOverride && !nextMaterialColorValid)
	this->markStale(STALE_SOURCE);
    if (changed) {
	(void)this->setCompactInstanceDisplayStateForPath("", TRUE,
	    visibilityChanged, nextVisible,
	    selectedChanged && !this->isCompactOccurrenceRegistry(),
	    nextSelected,
	    highlightedChanged && !this->isCompactOccurrenceRegistry(),
	    nextHighlighted);
	if (appearanceChanged && this->d->compactIndex) {
	    this->rebuildCompactInstanceDisplayState(TRUE);
	    this->markCompiledAssemblyDirty();
	    this->touch();
	}
	this->markCadBatchDirty();
	this->syncRealizedShapeOwnerState();
    }
    return changed;
}

int
SoBRLDatabaseSource::applyDisplayPatch(
    const BObolDatabaseSourceDisplayPatch &patch)
{
    int changed = 0;
    int appearanceChanged = 0;
    int visibilityChanged = 0;
    int selectedChanged = 0;
    int highlightedChanged = 0;
    if (patch.visibleValid && this->visible.getValue() != patch.visible) {
	this->visible = patch.visible;
	visibilityChanged = 1;
	changed = 1;
    }
    if (patch.selectedValid && this->selected.getValue() != patch.selected) {
	this->selected = patch.selected;
	selectedChanged = 1;
	changed = 1;
    }
    if (patch.highlightedValid &&
	this->highlighted.getValue() != patch.highlighted) {
	this->highlighted = patch.highlighted;
	highlightedChanged = 1;
	changed = 1;
    }
    if (patch.lineStyleValid &&
	this->lineStyle.getValue() != patch.lineStyle) {
	this->lineStyle = patch.lineStyle;
	appearanceChanged = 1;
	changed = 1;
    }
    if (patch.lineWidthValid &&
	this->lineWidth.getValue() != patch.lineWidth) {
	this->lineWidth = patch.lineWidth;
	appearanceChanged = 1;
	changed = 1;
    }
    if (patch.transparencyValid &&
	database_source_float_different(this->transparency.getValue(),
					patch.transparency)) {
	this->transparency = patch.transparency;
	appearanceChanged = 1;
	changed = 1;
    }
    if (patch.colorOverrideValid &&
	this->colorOverride.getValue() != patch.colorOverride) {
	this->colorOverride = patch.colorOverride;
	appearanceChanged = 1;
	changed = 1;
    }
    if (patch.colorValid &&
	!database_source_color_equal(this->color.getValue(), patch.color)) {
	this->color = patch.color;
	appearanceChanged = 1;
	changed = 1;
    }
    if (patch.selectedColorValid &&
	!database_source_color_equal(this->selectedColor.getValue(),
				     patch.selectedColor)) {
	this->selectedColor = patch.selectedColor;
	appearanceChanged = 1;
	changed = 1;
    }
    if (patch.highlightedColorValid &&
	!database_source_color_equal(this->highlightedColor.getValue(),
				     patch.highlightedColor)) {
	this->highlightedColor = patch.highlightedColor;
	appearanceChanged = 1;
	changed = 1;
    }
    if (patch.ghostedColorValid &&
	!database_source_color_equal(this->ghostedColor.getValue(),
				     patch.ghostedColor)) {
	this->ghostedColor = patch.ghostedColor;
	appearanceChanged = 1;
	changed = 1;
    }
    if (changed) {
	(void)this->setCompactInstanceDisplayStateForPath("", TRUE,
	    visibilityChanged, patch.visible,
	    selectedChanged && !this->isCompactOccurrenceRegistry(),
	    patch.selected,
	    highlightedChanged && !this->isCompactOccurrenceRegistry(),
	    patch.highlighted);
	if (appearanceChanged && this->d->compactIndex) {
	    this->rebuildCompactInstanceDisplayState(TRUE);
	    this->markCompiledAssemblyDirty();
	    this->touch();
	}
	this->markCadBatchDirty();
	this->syncRealizedShapeOwnerState();
    }
    return changed;
}

int
SoBRLDatabaseSource::setPlacementState(SbBool nextDrawMatrixValid,
				       const SbMatrix &nextDrawMatrix,
				       SbBool nextDrawCenterValid,
				       const SbVec3f &nextDrawCenter,
				       SbBool nextDrawSizeValid,
				       float nextDrawSize)
{
    int changed = 0;
    if (this->drawMatrixValid.getValue() != nextDrawMatrixValid) {
	this->drawMatrixValid = nextDrawMatrixValid;
	changed = 1;
    }
    if (!this->drawMatrix.getValue().equals(nextDrawMatrix, 0.000001f)) {
	this->drawMatrix = nextDrawMatrix;
	changed = 1;
    }
    if (this->drawCenterValid.getValue() != nextDrawCenterValid) {
	this->drawCenterValid = nextDrawCenterValid;
	changed = 1;
    }
    const SbVec3f currentCenter = this->drawCenter.getValue();
    if (database_source_float_different(currentCenter[0],
					nextDrawCenter[0]) ||
	database_source_float_different(currentCenter[1],
					nextDrawCenter[1]) ||
	database_source_float_different(currentCenter[2],
					nextDrawCenter[2])) {
	this->drawCenter = nextDrawCenter;
	changed = 1;
    }
    if (this->drawSizeValid.getValue() != nextDrawSizeValid) {
	this->drawSizeValid = nextDrawSizeValid;
	changed = 1;
    }
    if (database_source_float_different(this->drawSize.getValue(),
					nextDrawSize)) {
	this->drawSize = nextDrawSize;
	changed = 1;
    }
    if (sync_source_placement_transform(this))
	changed = 1;
    if (changed) {
	this->syncCompactInstancePlacementState();
	sync_realized_shape_placement_state_in_node(this, this);
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
	this->markDisplayMeshLodDirty();
    }
    return changed;
}

int
SoBRLDatabaseSource::setSourceBoundsState(SbBool nextBoundsValid,
	const SbVec3f &nextBoundsMin,
	const SbVec3f &nextBoundsMax)
{
    SbVec3f sanitizedMin(0.0f, 0.0f, 0.0f);
    SbVec3f sanitizedMax(0.0f, 0.0f, 0.0f);
    if (nextBoundsValid) {
	const SbBox3f bounds =
	    database_source_box_from_minmax(nextBoundsMin, nextBoundsMax);
	sanitizedMin = bounds.getMin();
	sanitizedMax = bounds.getMax();
    }

    int changed = 0;
    if (this->sourceBoundsValid.getValue() != nextBoundsValid) {
	this->sourceBoundsValid = nextBoundsValid;
	changed = 1;
    }
    if (!database_source_vec3f_equal(this->sourceBoundsMin.getValue(),
				     sanitizedMin)) {
	this->sourceBoundsMin = sanitizedMin;
	changed = 1;
    }
    if (!database_source_vec3f_equal(this->sourceBoundsMax.getValue(),
				     sanitizedMax)) {
	this->sourceBoundsMax = sanitizedMax;
	changed = 1;
    }

    return changed;
}

void
SoBRLDatabaseSource::clearSourceBounds(void)
{
    (void)this->setSourceBoundsState(FALSE,
				     SbVec3f(0.0f, 0.0f, 0.0f),
				     SbVec3f(0.0f, 0.0f, 0.0f));
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

void
SoBRLDatabaseSource::setDatabase(struct db_i *database)
{
    if (this->d->dbip != database) {
	this->d->dbip = database;
	this->markStale(STALE_DATABASE);
    }
}

struct db_i *
SoBRLDatabaseSource::getDatabase(void) const {
    return this->d->dbip;
}

void
SoBRLDatabaseSource::setMeshLod(struct BObolMeshLod *lod)
{
    if (this->d->meshLod == lod) {
	if (!lod)
	    this->clearMeshLodBounds();
	return;
    }

    if (this->d->meshLod)
	bobol_mesh_lod_destroy(this->d->meshLod);

    this->d->meshLod = lod;
    this->clearMeshLodBounds();
}

struct BObolMeshLod *
SoBRLDatabaseSource::getMeshLod(void) const {
    return this->d->meshLod;
}

void
SoBRLDatabaseSource::clearMeshLod(void)
{
    if (this->d->meshLod)
	bobol_mesh_lod_destroy(this->d->meshLod);
    this->d->meshLod = NULL;
    this->clearMeshLodBounds();
}

int
SoBRLDatabaseSource::setMeshLodBounds(const SbVec3f &bmin,
				      const SbVec3f &bmax)
{
    if (bmin[0] > bmax[0] || bmin[1] > bmax[1] || bmin[2] > bmax[2])
	return 0;

    this->d->meshLodBoundsMin = bmin;
    this->d->meshLodBoundsMax = bmax;
    this->d->meshLodBoundsValid = TRUE;
    return 1;
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
SoBRLDatabaseSource::clearMeshLodBounds(void)
{
    this->d->meshLodBoundsValid = FALSE;
    this->d->meshLodBoundsMin.setValue(0.0f, 0.0f, 0.0f);
    this->d->meshLodBoundsMax.setValue(0.0f, 0.0f, 0.0f);
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
    int sanitizedMode = mode == SHADED ? SHADED : WIREFRAME;
    uint32_t reason = STALE_NONE;
    if (this->d->dbip != database)
	reason |= STALE_DATABASE;
    if (this->drawMode.getValue() != sanitizedMode)
	reason |= STALE_DRAW;
    const std::string stableSourcePath = sourcePath ? sourcePath : "";
    const std::string stableInstanceKey =
	(sourceInstanceKey && sourceInstanceKey[0]) ?
	sourceInstanceKey : stableSourcePath.c_str();
    const char *effectiveInstanceKey = stableInstanceKey.c_str();
    if (bu_strcmp(this->instanceKey.getValue().getString(),
	       effectiveInstanceKey) != 0)
	reason |= STALE_SOURCE;
    const std::string stableRepresentationKey =
	(sourceRepresentationKey && sourceRepresentationKey[0]) ?
	sourceRepresentationKey : stableInstanceKey.c_str();
    const char *effectiveRepresentationKey = stableRepresentationKey.c_str();
    if (bu_strcmp(this->representationKey.getValue().getString(),
	       effectiveRepresentationKey) != 0)
	reason |= STALE_SOURCE;
    if (this->representationMode.getValue() != sourceRepresentationMode)
	reason |= STALE_DRAW;
    if (bu_strcmp(this->path.getValue().getString(),
	       stableSourcePath.c_str()) != 0)
	reason |= STALE_SOURCE;
    if (this->sourceRevision.getValue() != revision)
	reason |= STALE_SOURCE;

    this->detachFieldSensors();
    this->d->dbip = database;
    this->instanceKey = effectiveInstanceKey;
    this->path = stableSourcePath.c_str();
    this->representationKey = effectiveRepresentationKey;
    this->representationMode = sourceRepresentationMode;
    this->auxiliarySource = FALSE;
    this->drawMode = sanitizedMode;
    this->sourceRevision = revision;
    this->markStale(reason);
    this->attachFieldSensors();
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

SoBRLDatabaseSource *
SoBRLDatabaseSource::createDetachedRealizationTemplate(void) const
{
    SoBRLDatabaseSource *detached = new SoBRLDatabaseSource;
    detached->ref();
    detached->detachFieldSensors();
    detached->copyFieldValues(this, FALSE);
    detached->d->dbip = NULL;
    detached->stale = TRUE;
    detached->staleReason = STALE_SOURCE;
    detached->realizationStatus = UNREALIZED;
    detached->realizationDiagnostic = "";
    detached->attachFieldSensors();
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
     * A file-backed database is already the serialized input needed by
     * realization.  Give the worker its own read-only directory/file handle
     * instead of copying a multi-gigabyte draw closure into an in-memory DB
     * before it may publish one box.  Scene/source revisions still gate every
     * streamed merge and final adoption, so an edit that races this reader
     * invalidates its result.  Independent FILE state also avoids sharing the
     * live handle's seek cursor.
     *
     * In-memory databases and unusual file-open failures retain the closure
     * snapshot fallback: there is no separately reopenable serialized source
     * in those cases.
     */
    struct db_i *database = NULL;
    if (sourceDatabase->dbi_filename &&
	sourceDatabase->dbi_filename[0]) {
	database = db_open(sourceDatabase->dbi_filename, DB_OPEN_READONLY);
	if (database && db_dirbuild(database) < 0) {
	    db_close(database);
	    database = NULL;
	}
	if (database) {
	    (void)db_comb_instance_ids_set(database, 1);
	    database->dbi_read_only = 1;
	}
    }
    if (!database)
	database = database_snapshot_create(sourceDatabase,
	    this->path.getValue(), snapshotPathOut);
    if (!database)
	return FALSE;
    this->d->dbip = database;
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
    struct db_i *database = database_snapshot_create(this->d->dbip,
	this->path.getValue(), snapshotPathOut);
    if (!database) {
	detached->unref();
	return NULL;
    }
    detached->d->dbip = database;
    *databaseOut = database;
    return detached;
}

int
SoBRLDatabaseSource::adoptDetachedCompactRealization(
    SoBRLDatabaseSource *detached,
    SbBool authoritativeStreamDrained)
{
    if (!detached ||
	bu_strcmp(this->instanceKey.getValue().getString(),
	    detached->instanceKey.getValue().getString()) != 0 ||
	bu_strcmp(this->path.getValue().getString(),
	    detached->path.getValue().getString()) != 0 ||
	this->sourceRevision.getValue() !=
	    detached->sourceRevision.getValue() ||
	this->inputsRevision.getValue() !=
	    detached->inputsRevision.getValue() ||
	this->drawMode.getValue() != detached->drawMode.getValue() ||
	this->representationMode.getValue() !=
	    detached->representationMode.getValue())
	return 0;

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
    bool streamedComplete = authoritative ?
	(streamedAuthoritativeCount == authoritativeCount) :
	(authoritativeStreamDrained && streamedAuthoritativeCount > 0);
    if (streamedComplete && !authoritativeStreamDrained) {
	for (const BObolCompactInstanceEntry &entry :
	     authoritative->entries) {
	    const char *entryPath = database_source_skip_leading_slash(
		entry.semantic.path.getString());
	    const auto found = entryPath && entryPath[0] ?
		current->entryIndexByPath.find(entryPath) :
		current->entryIndexByPath.end();
	    if (found == current->entryIndexByPath.end() ||
		found->second >= current->entries.size()) {
		streamedComplete = false;
		break;
	    }
	    const BObolCompactInstanceEntry &live =
		current->entries[found->second];
	    if (live.sourceMeshRequestValid !=
		    entry.sourceMeshRequestValid ||
		bu_strcmp(live.shapeSummary.geometryKind.getString(),
		    entry.shapeSummary.geometryKind.getString()) != 0 ||
		!live.localToSource.equals(entry.localToSource, 0.000001f)) {
		streamedComplete = false;
		break;
	    }
	}
    }
    if (streamedComplete) {
	std::vector<size_t> retiredOverviewEntries;
	for (size_t i = 0; i < current->entries.size(); ++i) {
	    BObolCompactInstanceEntry &entry = current->entries[i];
	    if (!BU_STR_EQUAL(entry.shapeSummary.recordRole.getString(),
		    "lod-overview"))
		continue;
	    /*
	     * This is a publication barrier, not merely an index-state edit.
	     * Visibility-frontier maintenance may already have marked the
	     * synthetic overview invisible in the compact index while an
	     * earlier retained presentation is still displaying it.  Skipping
	     * an already-invisible entry therefore strands the box forever.
	     *
	     * Advance the revision unconditionally and journal this one record.
	     * Each view will then patch its hidden set even when the index value
	     * itself was already false.  Keeping the terminal delta sparse also
	     * avoids turning retirement of one tiny startup aid into a 50k
	     * occurrence presentation scan/rebuild.
	     */
	    entry.authoredVisible = FALSE;
	    entry.visible = FALSE;
	    entry.visibilityRevision =
		compact_next_revision(entry.visibilityRevision);
	    compact_sync_shape_summary_state(entry);
	    current->hiddenInstances.push_back(entry.instance);
	    retiredOverviewEntries.push_back(i);
	}
	if (!retiredOverviewEntries.empty()) {
	    std::sort(current->hiddenInstances.begin(),
		current->hiddenInstances.end());
	    current->hiddenInstances.erase(std::unique(
		current->hiddenInstances.begin(),
		current->hiddenInstances.end()),
		current->hiddenInstances.end());
	    this->markCompiledAssemblyDirty();
	    this->markCadBatchDirty(retiredOverviewEntries);
	}
	remove_non_auxiliary_children(this);

	SbBox3f bounds;
	if ((detached->getSourceBounds(bounds) ||
		this->getSourceBounds(bounds)) && !bounds.isEmpty())
	    (void)this->setSourceBoundsState(TRUE, bounds.getMin(),
		bounds.getMax());
	else
	    this->clearSourceBounds();

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

	this->realizedRevision = this->sourceRevision.getValue();
	this->realizedSourceRevision = this->sourceRevision.getValue();
	this->realizedInputsRevision = this->inputsRevision.getValue();
	this->realizedViewRevision = this->viewRevision.getValue();
	this->realizationStatus = REALIZED;
	this->realizationDiagnostic = "";
	this->realizationIdentity = source_realization_identity(this);
	this->stale = FALSE;
	this->staleReason = STALE_NONE;
	if (getenv("BOBOL_DRAW_TIMING"))
	    bu_log("[obol-timing] deferred adoption preserved streamed "
		   "index: n=%zu overview=%zu\n",
		   authoritativeCount, overviewCount);
	return authoritativeCount > static_cast<size_t>(INT_MAX) ? INT_MAX :
	    static_cast<int>(authoritativeCount);
    }

    if (!authoritative)
	return 0;

    /*
     * Camera revision is not source-content identity.  A user commonly moves
     * the view while this worker replaces a structural frontier with retained
     * native/PoP-backed occurrences.  Rejecting that valid result strands
     * depth-capped combination boxes; the active PoP prefix is selected from
     * the current view independently after adoption.
     */
    BObolCompactInstanceIndex *next = detached->d->compactIndex;
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
    this->d->compactHandleSourceId = detached->d->compactHandleSourceId;
    this->installCompactInstanceIndex(next, TRUE);
    remove_non_auxiliary_children(this);
    this->markCompiledAssemblyDirty();
    this->markCadBatchDirty();

    SbBox3f bounds;
    if (detached->getSourceBounds(bounds) && !bounds.isEmpty())
	(void)this->setSourceBoundsState(TRUE, bounds.getMin(), bounds.getMax());
    else
	this->clearSourceBounds();

    this->realizedRevision = this->sourceRevision.getValue();
    this->realizedSourceRevision = this->sourceRevision.getValue();
    this->realizedInputsRevision = this->inputsRevision.getValue();
    this->realizedViewRevision = this->viewRevision.getValue();
    this->realizationStatus = REALIZED;
    this->realizationDiagnostic = "";
    this->realizationIdentity = source_realization_identity(this);
    this->stale = FALSE;
    this->staleReason = STALE_NONE;
    return static_cast<int>(this->d->compactIndex->entries.size());
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

    database_source_assign_string(shape->sourcePath, newSourcePath);
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

    SbString identity = source_record_identity(source, newSourcePath);
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
    object->sourcePath = newSourcePath;
    object->sourceId = revision;
    if (!sourceName || !sourceName[0] ||
	(oldName && BU_STR_EQUAL(sourceName, oldName)) ||
	(oldSourcePath && BU_STR_EQUAL(sourceName, oldSourcePath)))
	object->sourceName = newName;
}

static void
retarget_realized_node_source(SoNode *node,
			      const SoBRLDatabaseSource *source,
			      const char *oldSourcePath,
			      const char *newSourcePath,
			      uint32_t revision)
{
    if (!node || !newSourcePath)
	return;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	retarget_realized_shape_source(static_cast<SoBRLVListShape *>(node),
				       source, oldSourcePath, newSourcePath, revision);
	return;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	retarget_realized_shape_source(static_cast<SoBRLMeshShape *>(node),
				       source, oldSourcePath, newSourcePath, revision);
	return;
    }
    if (node->isOfType(SoBRLMaterialObject::getClassTypeId())) {
	retarget_material_object_source(
	    static_cast<SoBRLMaterialObject *>(node), oldSourcePath,
	    newSourcePath, revision);
	return;
    }
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    retarget_realized_node_source(group->getChild(i), source,
					  oldSourcePath, newSourcePath, revision);
    }
}

int
SoBRLDatabaseSource::retargetDatabaseSource(const char *sourcePath,
	uint32_t revision)
{
    return this->retargetDatabaseSourceInstance(sourcePath, sourcePath,
	    revision);
}

int
SoBRLDatabaseSource::retargetDatabaseSourceInstance(
    const char *sourceInstanceKey,
    const char *sourcePath,
    uint32_t revision)
{
    if (!sourcePath || !sourcePath[0])
	return -1;

    const std::string stableSourcePath = sourcePath;
    const std::string stableInstanceKey =
	(sourceInstanceKey && sourceInstanceKey[0]) ?
	sourceInstanceKey : stableSourcePath.c_str();
    const char *effectiveInstanceKey = stableInstanceKey.c_str();
    const char *oldPath = this->path.getValue().getString();
    const char *oldInstanceKey = this->instanceKey.getValue().getString();
    const int pathChanged = !oldPath ||
			    !BU_STR_EQUAL(oldPath, stableSourcePath.c_str());
    const int instanceChanged = !oldInstanceKey ||
				!BU_STR_EQUAL(oldInstanceKey, effectiveInstanceKey);
    if (revision == 0)
	revision = this->sourceRevision.getValue() +
		   ((pathChanged || instanceChanged) ? 1 : 0);
    const int revisionChanged = this->sourceRevision.getValue() != revision;
    if (!pathChanged && !instanceChanged && !revisionChanged)
	return 0;
    const int sourceChanged = pathChanged || revisionChanged;

    std::string oldPathCopy = oldPath ? oldPath : "";

    this->detachFieldSensors();
    this->instanceKey = effectiveInstanceKey;
    this->path = stableSourcePath.c_str();
    this->sourceRevision = revision;
    if (sourceChanged)
	this->markStale(STALE_SOURCE);
    for (int i = 0; i < this->getNumChildren(); i++)
	retarget_realized_node_source(this->getChild(i), this,
				      oldPathCopy.c_str(), stableSourcePath.c_str(), revision);
    this->realizationIdentity = source_realization_identity(this);
    this->attachFieldSensors();
    return 1;
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
bobol_database_source_realize_wireframe_compact_with_cache(
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
	    revision);
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
	mark_source_realized_current(source);
	return 1;
    }
    if (directRealized < 0) {
	if (!cache->preserveCompactSourceOnFailure) {
	    remove_non_auxiliary_children(source);
	    source->discardCompactInstanceHistory();
	    source->realizationIdentity = "";
	}
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
    data.compact_index = new BObolCompactInstanceIndex;
    data.stream_sink = stream;
    data.material_sweep = &materialSweep;

    const char *av[1] = {treeName};
    const int ret = db_walk_tree_leaf_instances(source->d->dbip, 1, av, 1,
	&init_state, NULL, NULL, realize_leaf, &data);
    db_free_db_tree_state(&init_state);

    if (ret < 0 || data.realized_shapes <= 0 || data.failed_shapes > 0 ||
	data.compact_unsupported || data.compact_index->entries.empty()) {
	const size_t compactEntryCount = data.compact_index->entries.size();
	delete data.compact_index;
	if (!cache->preserveCompactSourceOnFailure) {
	    remove_non_auxiliary_children(source);
	    source->discardCompactInstanceHistory();
	    source->realizationIdentity = "";
	}
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

    source->installCompactInstanceIndex(data.compact_index, TRUE);
    source->markCompiledAssemblyDirty();
    if (data.compact_bounds_valid && !data.compact_bounds.isEmpty())
	(void)source->setSourceBoundsState(TRUE, data.compact_bounds.getMin(),
	    data.compact_bounds.getMax());
    else
	source->clearSourceBounds();
    mark_source_realized_current(source);
    bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_SOURCES, 1);
    bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_INSTANCES,
	static_cast<uint64_t>(source->d->compactIndex->entries.size()));
    return 1;
}

SbBool
bobol_database_source_realize_wireframe_with_cache(
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
	    update_source_bounds_from_realized_geometry(source);
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
	source->realizedRevision = source->sourceRevision.getValue();
	source->realizedSourceRevision = source->sourceRevision.getValue();
	source->realizedInputsRevision = source->inputsRevision.getValue();
	source->realizedViewRevision = source->viewRevision.getValue();
	source->realizationStatus = SoBRLDatabaseSource::REALIZED;
	source->realizationDiagnostic = "";
	source->realizationIdentity = source_realization_identity(source);
	source->stale = FALSE;
	source->staleReason = SoBRLDatabaseSource::STALE_NONE;
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

    source->realizedRevision = source->sourceRevision.getValue();
    source->realizedSourceRevision = source->sourceRevision.getValue();
    source->realizedInputsRevision = source->inputsRevision.getValue();
    source->realizedViewRevision = source->viewRevision.getValue();
    source->realizationStatus = SoBRLDatabaseSource::REALIZED;
    source->realizationDiagnostic = "";
    source->realizationIdentity = source_realization_identity(source);
    source->stale = FALSE;
    source->staleReason = SoBRLDatabaseSource::STALE_NONE;
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
    bool coverageReady = false;
    /* True only when geometry is the standing proxy for a source-backed
     * triangle PoP asset.  BREP wire drawing has its own immutable
     * progressive line representation and must not also submit the shaded
     * tessellation as a wire overlay. */
    bool sourceMeshReady = false;
    bool ready = false;
    std::once_flag coverageOnce;
    std::once_flag realizeOnce;
};

struct compact_coverage_work_item {
    compact_coverage_asset *asset = NULL;
    compact_coverage_occurrence occurrence;
};

struct compact_coverage_collect {
    SoBRLDatabaseSource *source = NULL;
    BObolDatabaseSourceRealizationCache *cache = NULL;
    BObolCompactOccurrenceStream *stream = NULL;
    BObolMaterialColorSweep *materialSweep = NULL;
    uint32_t revision = 0;
    std::vector<std::unique_ptr<compact_coverage_asset>> assets;
    std::unordered_map<std::string, size_t> assetIndices;
    std::unordered_set<std::string> seenInstances;
    std::unordered_map<std::string, uint32_t> occurrenceCounts;
    size_t occurrenceCount = 0;
    std::mutex workMutex;
    std::condition_variable workReady;
    std::deque<compact_coverage_work_item> work;
    /* Coverage is the latency-critical phase.  Retain its completed work
     * records here, then import full BoTs only after every leaf box and the
     * exact target extent have been published. */
    std::deque<compact_coverage_work_item> detailWork;
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
 * db_get_external still owns one encoded-object copy.  The independent
 * read-only realization database is memory mapped, making these calls
 * parallel and free of the live editor handle's shared stdio cursor.  V4
 * retains the complete-import fallback because its record layout does not
 * offer the same compact body contract.
 */
static bool
compact_coverage_bot_bounds(struct db_i *dbip, struct directory *dp,
	SbBox3f &bounds, size_t &vertexCount, size_t &faceCount,
	unsigned char &mode, unsigned char &orientation)
{
    bounds.makeEmpty();
    vertexCount = 0;
    faceCount = 0;
    mode = 0;
    orientation = 0;
    if (!dbip || !dp ||
	dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT)
	return false;

    if (db_version(dbip) != 5) {
	struct rt_db_internal intern;
	RT_DB_INTERNAL_INIT(&intern);
	if (rt_db_get_internal(&intern, dp, dbip, NULL) < 0)
	    return false;
	const struct rt_bot_internal *bot =
	    intern.idb_type == ID_BOT && intern.idb_ptr ?
	    static_cast<const struct rt_bot_internal *>(intern.idb_ptr) :
	    NULL;
	if (bot) {
	    BObolSourceMeshRequest request;
	    if (cad_source_mesh_request_from_bot(request, bot)) {
		bounds = request.bounds;
		vertexCount = bot->num_vertices;
		faceCount = bot->num_faces;
		mode = bot->mode;
		orientation = bot->orientation;
	    }
	}
	rt_db_free_internal(&intern);
	return !bounds.isEmpty();
    }

    struct bu_external external;
    BU_EXTERNAL_INIT(&external);
    if (db_get_external(&external, dp, dbip) < 0)
	return false;

    struct db5_raw_internal raw;
    const bool rawValid =
	db5_get_raw_internal_ptr(&raw, external.ext_buf) != NULL &&
	raw.major_type == DB5_MAJORTYPE_BRLCAD &&
	raw.minor_type == DB5_MINORTYPE_BRLCAD_BOT &&
	raw.body.ext_buf && raw.body.ext_nbytes >=
	    2 * SIZEOF_NETWORK_LONG + 3;
    if (!rawValid) {
	bu_free_external(&external);
	return false;
    }

    const unsigned char *body = raw.body.ext_buf;
    const size_t vertices = static_cast<size_t>(BU_GLONG(body));
    const size_t faces = static_cast<size_t>(
	BU_GLONG(body + SIZEOF_NETWORK_LONG));
    const size_t fixedBytes = 2 * SIZEOF_NETWORK_LONG + 3;
    const size_t vertexStride =
	SIZEOF_NETWORK_DOUBLE * ELEMENTS_PER_POINT;
    if (!vertices || vertices > (SIZE_MAX - fixedBytes) / vertexStride ||
	fixedBytes + vertices * vertexStride > raw.body.ext_nbytes) {
	bu_free_external(&external);
	return false;
    }

    orientation = body[2 * SIZEOF_NETWORK_LONG];
    mode = body[2 * SIZEOF_NETWORK_LONG + 1];
    const unsigned char *point = body + fixedBytes;
    point_t minimum = {INFINITY, INFINITY, INFINITY};
    point_t maximum = {-INFINITY, -INFINITY, -INFINITY};
    size_t finiteVertices = 0;
    for (size_t i = 0; i < vertices; ++i) {
	double value[ELEMENTS_PER_POINT] = {0.0, 0.0, 0.0};
	bu_cv_ntohd(reinterpret_cast<unsigned char *>(value), point,
	    ELEMENTS_PER_POINT);
	point += vertexStride;
	if (!std::isfinite(value[X]) || !std::isfinite(value[Y]) ||
	    !std::isfinite(value[Z]))
	    continue;
	for (size_t axis = 0; axis < 3; ++axis) {
	    minimum[axis] = std::min(minimum[axis], value[axis]);
	    maximum[axis] = std::max(maximum[axis], value[axis]);
	}
	finiteVertices++;
    }
    bu_free_external(&external);
    if (!finiteVertices)
	return false;

    bounds = SbBox3f(
	SbVec3f(static_cast<float>(minimum[X]),
	    static_cast<float>(minimum[Y]),
	    static_cast<float>(minimum[Z])),
	SbVec3f(static_cast<float>(maximum[X]),
	    static_cast<float>(maximum[Y]),
	    static_cast<float>(maximum[Z])));
    vertexCount = vertices;
    faceCount = faces;
    return !bounds.isEmpty();
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
    if (dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT &&
	dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BREP)
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
    auto foundAsset = collect->assetIndices.find(cacheKey);
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
	collect->assetIndices[cacheKey] = assetIndex;
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
	dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP ? "brep" : "bot";
    occurrence.summary = compact_occurrence_tree_summary(
	collect->source, tsp, semanticPath.c_str(), dp->d_namep, sourceType,
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
    compact_coverage_work_item work;
    work.asset = collect->assets[assetIndex].get();
    work.occurrence = std::move(occurrence);
    {
	std::lock_guard<std::mutex> guard(collect->workMutex);
	collect->work.push_back(std::move(work));
	collect->occurrenceCount++;
    }
    collect->workReady.notify_one();
    bu_free(rawPath, "compact coverage path");
    return make_nop_tree();
}

static std::string
compact_coverage_broad_key(const compact_coverage_asset &asset)
{
    if (!asset.ready || !asset.sampleFingerprintValid)
	return std::string();
    char key[192] = {0};
    snprintf(key, sizeof(key), "%zu:%zu:%u:%u:%016llx",
	asset.vertexCount, asset.faceCount,
	static_cast<unsigned int>(asset.mode),
	static_cast<unsigned int>(asset.orientation),
	static_cast<unsigned long long>(asset.sampleFingerprint));
    return std::string(key);
}

static std::shared_ptr<const Obol::PartGeometry>
compact_coverage_aabb_geometry(const SbBox3f &bounds,
	SbMatrix &geometryTransform)
{
    geometryTransform = SbMatrix::identity();
    if (bounds.isEmpty())
	return std::shared_ptr<const Obol::PartGeometry>();
    const SbVec3f minimum = bounds.getMin();
    const SbVec3f maximum = bounds.getMax();
    const SbVec3f extent = maximum - minimum;
    if (extent[0] > SMALL_FASTF && extent[1] > SMALL_FASTF &&
	extent[2] > SMALL_FASTF) {
	static const std::shared_ptr<const Obol::PartGeometry> unitAabb = []() {
	    Obol::PartGeometry geometry;
	    const SbBox3f unitBounds(SbVec3f(0.0f, 0.0f, 0.0f),
		SbVec3f(1.0f, 1.0f, 1.0f));
	    if (!cad_wire_part_geometry_from_aabb(unitBounds, geometry))
		return std::shared_ptr<const Obol::PartGeometry>();
	    geometry.subpixelProxyEligible = true;
	    geometry.structuralProxy = true;
	    return std::make_shared<const Obol::PartGeometry>(
		std::move(geometry));
	}();
	geometryTransform.setScale(extent);
	SbMatrix translation;
	translation.setTranslate(minimum);
	geometryTransform.multRight(translation);
	return unitAabb;
    }
    Obol::PartGeometry geometry;
    if (!cad_wire_part_geometry_from_aabb(bounds, geometry))
	return std::shared_ptr<const Obol::PartGeometry>();
    return std::make_shared<const Obol::PartGeometry>(
	std::move(geometry));
}

static BObolCompactOccurrence
compact_coverage_overview_occurrence(SoBRLDatabaseSource *source,
	const char *treeName, const SbBox3f &bounds, uint32_t revision)
{
    BObolCompactOccurrence overview;
    if (!source || !treeName || !treeName[0] || bounds.isEmpty())
	return overview;

    overview.geometry = compact_coverage_aabb_geometry(bounds,
	overview.geometryTransform);
    if (!overview.geometry)
	return overview;
    std::string overviewPath(treeName);
    overviewPath += "/@draw-extent";
    overview.summary = compact_occurrence_summary(source,
	overviewPath.c_str(), treeName, "proxy", "overview-aabb", revision,
	BObolRealizedShapeSummary::SHAPE_VLIST);
    overview.summary.recordRole = "lod-overview";
    overview.summary.selectable = FALSE;
    overview.summary.lodAvailable = TRUE;
    overview.summary.lodActiveLevel = BOBOL_LOD_QUALITY_PROXY;
    overview.summary.lodBoundsMin = bounds.getMin();
    overview.summary.lodBoundsMax = bounds.getMax();
    overview.summary.pointCount = 24;
    overview.summary.commandCount = 24;
    overview.summary.segmentCount = 12;
    overview.summary.boundsValid = TRUE;
    overview.summary.bounds = bounds;
    return overview;
}

static int
compact_stream_publish_parallel_coverage(
	SoBRLDatabaseSource *source,
	BObolDatabaseSourceRealizationCache *cache,
	const char *treeName,
	BObolCompactOccurrenceStream *stream,
	uint32_t revision)
{
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
    std::atomic<size_t> publishedBoxes(0);
    std::atomic<size_t> publishedContracts(0);
    std::mutex aggregateMutex;
    SbBox3f aggregateBounds;
    aggregateBounds.makeEmpty();
    int64_t lastOverviewPublication = 0;
    size_t workerCount = bu_avail_cpus();
    workerCount = std::max<size_t>(1, std::min<size_t>(workerCount, 32));
    const auto coverageWorker = [&]() {
	bu_nice_set(5);
	for (;;) {
	    compact_coverage_work_item item;
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
		item = std::move(collect.work.front());
		collect.work.pop_front();
	    }
	    if (!item.asset)
		continue;
	    compact_coverage_asset &asset = *item.asset;
	    std::call_once(asset.coverageOnce, [&]() {
		size_t coverageWorkingSet = asset.estimatedWorkingSetBytes;
		if (db_version(source->getDatabase()) == 5 && asset.dp &&
		    asset.dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) {
		    const size_t fixedBytes = 1024ULL * 1024ULL;
		    const size_t encodedBytes = asset.dp->d_len;
		    coverageWorkingSet =
			encodedBytes > SIZE_MAX - fixedBytes ?
			SIZE_MAX : encodedBytes + fixedBytes;
		}
		bobol_lod_working_set_acquire(coverageWorkingSet);
		if (!stream->isCancelled() && asset.dp &&
		    asset.dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) {
		    asset.coverageReady = compact_coverage_bot_bounds(
			source->getDatabase(), asset.dp, asset.coverageBounds,
			asset.vertexCount, asset.faceCount, asset.mode,
			asset.orientation) &&
			asset.faceCount >= source->lodBotThreshold.getValue();
		} else if (!stream->isCancelled() && asset.dp &&
		    asset.dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP) {
		    asset.coverageReady = compact_coverage_primitive_bounds(
			source->getDatabase(), asset.dp, asset.coverageBounds);
		}
		if (asset.coverageReady) {
		    asset.coverageGeometry = compact_coverage_aabb_geometry(
			asset.coverageBounds, asset.proxyGeometryTransform);
		    asset.coverageReady =
			asset.coverageGeometry ? true : false;
		}
		bobol_lod_working_set_release(coverageWorkingSet);
	    });

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
		if (!overviewSnapshot.isEmpty()) {
		    BObolCompactOccurrence overview =
			compact_coverage_overview_occurrence(source, treeName,
			    overviewSnapshot, revision);
		    if (overview.geometry)
			stream->pushPriority(overview);
		}

		BObolCompactOccurrence occurrence;
		occurrence.geometry = asset.coverageGeometry;
		occurrence.geometryTransform =
		    asset.proxyGeometryTransform;
		occurrence.localTransform = item.occurrence.localTransform;
		occurrence.lodBacked = TRUE;
		occurrence.sourceMeshRequestValid = FALSE;
		occurrence.occurrenceIndex =
		    item.occurrence.occurrenceIndex;
		occurrence.booleanOperation =
		    item.occurrence.booleanOperation;
		occurrence.summary = item.occurrence.summary;
		occurrence.summary.lodAvailable = TRUE;
		occurrence.summary.lodActiveLevel =
		    BOBOL_LOD_QUALITY_PROXY;
		occurrence.summary.lodFaceCount = asset.faceCount;
		occurrence.summary.lodPointCount = asset.vertexCount;
		occurrence.summary.lodOriginalPointCount =
		    asset.vertexCount;
		occurrence.summary.lodBoundsMin =
		    asset.coverageBounds.getMin();
		occurrence.summary.lodBoundsMax =
		    asset.coverageBounds.getMax();
		occurrence.summary.boundsValid = TRUE;
		occurrence.summary.bounds = asset.coverageBounds;
		stream->push(std::move(occurrence));
		publishedBoxes.fetch_add(1);
	    }

	    /* Full source import is deliberately a second phase. */
	    {
		std::lock_guard<std::mutex> guard(collect.workMutex);
		collect.detailWork.push_back(std::move(item));
	    }
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
	    if (!item.asset || stream->isCancelled())
		continue;
	    compact_coverage_asset &asset = *item.asset;
	    std::call_once(asset.realizeOnce, [&]() {
		bobol_lod_working_set_acquire(
		    asset.estimatedWorkingSetBytes);
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
			    asset.geometry = compact_coverage_aabb_geometry(
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
		Obol::PartGeometry wireGeometry;
		if (cad_progressive_wire_part_geometry_from_provider(
			&intern, source, wireGeometry) > 0 &&
		    wireGeometry.wire) {
		    asset.vertexCount = wireGeometry.wire->segmentCount() * 2;
		    asset.faceCount = 0;
		    asset.proxyGeometryTransform = SbMatrix::identity();
		    asset.geometry =
			std::make_shared<const Obol::PartGeometry>(
			    std::move(wireGeometry));
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
			    compact_coverage_aabb_geometry(
				asset.sourceMeshRequest.bounds,
				asset.proxyGeometryTransform);
			asset.ready = asset.geometry ? true : false;
			asset.sourceMeshReady = asset.ready;
		    }
		    if (asset.ready && !stream->isCancelled() &&
			stream->retainStagedSource(staged))
			asset.sourceMeshRequest.stagedSource = staged;
		}
	    }
		}
		if (ownsInternal)
		    rt_db_free_internal(&intern);
		bobol_lod_working_set_release(
		    asset.estimatedWorkingSetBytes);
	    });

	    if (!asset.ready || stream->isCancelled())
		continue;

	    BObolCompactOccurrence occurrence;
	    occurrence.geometry = asset.geometry;
	    occurrence.geometryTransform =
		asset.proxyGeometryTransform;
	    occurrence.localTransform = item.occurrence.localTransform;
	    occurrence.lodBacked = TRUE;
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
		occurrence.summary.geometryKind = "wire";
		occurrence.summary.pointCount = asset.vertexCount;
		occurrence.summary.commandCount = asset.vertexCount;
		occurrence.summary.segmentCount =
		    asset.geometry->wire->segmentCount();
		occurrence.summary.boundsValid = TRUE;
		occurrence.summary.bounds = asset.geometry->wire->bounds;
	    }
	    stream->push(std::move(occurrence));
	    publishedContracts.fetch_add(1);
	}
    };
    std::vector<std::thread> workers;
    workers.reserve(workerCount);
    for (size_t i = 0; i < workerCount; i++)
	workers.push_back(std::thread(coverageWorker));

    /*
     * Enumeration is the producer, not a prerequisite.  The old two-phase
     * walk held every occurrence until the complete hierarchy was known,
     * creating seconds of blank UI on 50k+ cold scenes while all CPUs sat
     * idle.  Workers now import bounds and publish the expanding overview/
     * leaf boxes while db_walk_tree continues discovering later branches.
     */
    struct db_tree_state initialState;
    db_init_db_tree_state(&initialState, source->getDatabase());
    initialState.ts_stop_at_regions = 0;
    const char *treeNames[1] = {treeName};
    const int walkResult = db_walk_tree_leaf_instances(
	source->getDatabase(), 1, treeNames, 1, &initialState, NULL, NULL,
	compact_coverage_collect_leaf, &collect);
    db_free_db_tree_state(&initialState);
    /* Publish the final occurrence cardinality as soon as enumeration ends,
     * while bound workers are still draining.  The GUI can then reserve its
     * relocation-sensitive compact vectors once instead of repeatedly
     * growing them during the long tail of a 50k stream. */
    if (walkResult >= 0)
	stream->setExpectedCount(collect.occurrenceCount);
    {
	std::lock_guard<std::mutex> guard(collect.workMutex);
	collect.producerDone = true;
    }
    collect.workReady.notify_all();
    for (std::thread &thread : workers)
	thread.join();
    if (walkResult < 0 || stream->isCancelled())
	return -1;
    if (collect.assets.empty())
	return 0;

    /*
     * Leaf AABBs may already fill the producer queue faster than the GUI can
     * adopt them.  Once the parallel bound pass knows every occurrence,
     * publish one unselectable whole-target extent through the priority lane.
     * It is intentionally absent from the final authoritative compact index;
     * adoption therefore removes it atomically after all leaf boxes/meshes
     * are present.
     */
    if (!aggregateBounds.isEmpty()) {
	BObolCompactOccurrence overview =
	    compact_coverage_overview_occurrence(source, treeName,
		aggregateBounds, revision);
	if (overview.geometry) {
	    stream->pushPriority(overview);
	    /*
	     * The priority lane is drained before ordinary leaves.  Publishing
	     * this release flag after enqueue therefore guarantees that the next
	     * provider tick which observes it has already merged the exact
	     * extent and may perform one stable deferred autoview.
	     */
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

    /*
     * Bounds-first is a scheduling barrier, not a loss of parallelism.  Once
     * complete visual coverage and an exact target extent are available,
     * reuse the same bounded worker count to import the source arrays needed
     * by PoP and enrich each standing leaf in place.
     */
    workers.clear();
    for (size_t i = 0; i < workerCount; i++)
	workers.push_back(std::thread(detailWorker));
    for (std::thread &thread : workers)
	thread.join();
    if (stream->isCancelled())
	return -1;
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
	       publishedBoxes.load(), publishedContracts.load(), seeded,
	       bobol_lod_working_set_global_peak_bytes(),
	       bobol_lod_working_set_global_peak_tasks());
    return publishedBoxes.load() > static_cast<size_t>(INT_MAX) ? INT_MAX :
	static_cast<int>(publishedBoxes.load());
}

} /* anonymous namespace */

int
bobol_database_source_realize_mesh_compact_with_cache(
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
    int directRealized = 0;
    {
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
	mark_source_realized_current(source);
	return 1;
    }
    if (directRealized < 0) {
	if (!cache->preserveCompactSourceOnFailure) {
	    remove_non_auxiliary_children(source);
	    source->discardCompactInstanceHistory();
	    source->realizationIdentity = "";
	}
	return -1;
    }

    /* Prefill is a batch optimization for synchronous static realization.
     * It is the wrong latency and memory tradeoff for progressive BoTs: the
     * occurrence walk already performs bounded online PCA while structural
     * boxes remain visible.  Running both paths can import two candidate/
     * representative pairs concurrently. */
    if (!stream && source->lodBotThreshold.getValue() == 0 &&
	compact_mesh_prefill_cache(source, cache, treeName, NULL) < 0) {
	if (!cache->preserveCompactSourceOnFailure) {
	    remove_non_auxiliary_children(source);
	    source->discardCompactInstanceHistory();
	    source->realizationIdentity = "";
	}
	return -1;
    }

    if (stream && source->lodBotThreshold.getValue() > 0 &&
	!stream->hasWarmCoverageComplete() &&
	compact_stream_publish_parallel_coverage(source, cache, treeName,
	    stream, revision) < 0) {
	if (!cache->preserveCompactSourceOnFailure) {
	    remove_non_auxiliary_children(source);
	    source->discardCompactInstanceHistory();
	    source->realizationIdentity = "";
	}
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
    data.compact_index = new BObolCompactInstanceIndex;
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
	delete data.compact_index;
	if (!cache->preserveCompactSourceOnFailure) {
	    remove_non_auxiliary_children(source);
	    source->discardCompactInstanceHistory();
	    source->realizationIdentity = "";
	}
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

    source->installCompactInstanceIndex(data.compact_index, TRUE);
    source->markCompiledAssemblyDirty();
    if (data.compact_bounds_valid && !data.compact_bounds.isEmpty())
	(void)source->setSourceBoundsState(TRUE, data.compact_bounds.getMin(),
	    data.compact_bounds.getMax());
    else
	source->clearSourceBounds();
    mark_source_realized_current(source);
    bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_SOURCES, 1);
    bobol_performance_counter_add(BOBOL_PERF_CAD_COMPACT_INSTANCES,
	static_cast<uint64_t>(source->d->compactIndex->entries.size()));
    return 1;
}

SbBool
bobol_database_source_realize_mesh_with_cache(
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
	    update_source_bounds_from_realized_geometry(source);
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
	source->realizedRevision = source->sourceRevision.getValue();
	source->realizedSourceRevision = source->sourceRevision.getValue();
	source->realizedInputsRevision = source->inputsRevision.getValue();
	source->realizedViewRevision = source->viewRevision.getValue();
	source->realizationStatus = SoBRLDatabaseSource::REALIZED;
	source->realizationDiagnostic = "";
	source->realizationIdentity = source_realization_identity(source);
	source->stale = FALSE;
	source->staleReason = SoBRLDatabaseSource::STALE_NONE;
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

    source->realizedRevision = source->sourceRevision.getValue();
    source->realizedSourceRevision = source->sourceRevision.getValue();
    source->realizedInputsRevision = source->inputsRevision.getValue();
    source->realizedViewRevision = source->viewRevision.getValue();
    source->realizationStatus = SoBRLDatabaseSource::REALIZED;
    source->realizationDiagnostic = "";
    source->realizationIdentity = source_realization_identity(source);
    source->stale = FALSE;
    source->staleReason = SoBRLDatabaseSource::STALE_NONE;
    update_source_bounds_from_realized_geometry(source);
    source->syncRealizedShapeOwnerState();
    return TRUE;
}

SbBool
SoBRLDatabaseSource::realizePrototypeWireframe(void)
{
    this->realizationDiagnostic = "";
    SoBRLVListShape *shape = this->getRealizedShape();
    if (!shape) {
	shape = new SoBRLVListShape;
	database_source_add_realized_child(this, shape);
    }

    const float halfExtent = 1.0f + 0.25f * static_cast<float>(this->sourceRevision.getValue() % 4);
    SbVec3f points[5] = {
	SbVec3f(-halfExtent, -halfExtent, 0.0f),
	SbVec3f(halfExtent, -halfExtent, 0.0f),
	SbVec3f(halfExtent,  halfExtent, 0.0f),
	SbVec3f(-halfExtent,  halfExtent, 0.0f),
	SbVec3f(-halfExtent, -halfExtent, 0.0f)
    };
    int32_t commands[5] = {
	SoBRLVListShape::MOVE,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW
    };

    shape->sourcePath = this->path.getValue();
    shape->sourceName = lookup_name_from_path(this->path.getValue()) ?
			lookup_name_from_path(this->path.getValue()) : this->path.getValue().getString();
    shape->sourceType = "prototype";
    shape->sourceId = this->sourceRevision.getValue();
    shape->displayName = this->displayName.getValue().getLength() > 0 ?
			 this->displayName.getValue() : shape->sourceName.getValue();
    shape->geometryName = shape->sourceName.getValue();
    shape->sourceIdentity = source_record_identity(this,
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
    shape->setLineSet(points, commands, 5);

    this->realizedRevision = this->sourceRevision.getValue();
    this->realizedSourceRevision = this->sourceRevision.getValue();
    this->realizedInputsRevision = this->inputsRevision.getValue();
    this->realizedViewRevision = this->viewRevision.getValue();
    this->realizationStatus = REALIZED;
    this->realizationDiagnostic = "";
    this->realizationIdentity = source_realization_identity(this);
    this->stale = FALSE;
    this->staleReason = STALE_NONE;
    this->syncRealizedShapeOwnerState();
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

static void
set_external_bounds_from_points(SoBRLDatabaseSource *source,
				const SbVec3f *points,
				int count)
{
    if (!source)
	return;

    SbVec3f boundsMin;
    SbVec3f boundsMax;
    if (external_bounds_from_points(points, count, boundsMin, boundsMax))
	(void)source->setSourceBoundsState(TRUE, boundsMin, boundsMax);
    else
	source->clearSourceBounds();
}

static void
mark_external_primary_published_current(SoBRLDatabaseSource *source)
{
    if (!source)
	return;

    const uint32_t sourceRevision = source->sourceRevision.getValue();
    source->realizedRevision = sourceRevision;
    source->realizedSourceRevision = sourceRevision;
    source->realizedInputsRevision = source->inputsRevision.getValue();
    source->realizedViewRevision = source->viewRevision.getValue();
    source->realizationStatus = SoBRLDatabaseSource::REALIZED;
    source->realizationDiagnostic = "";
    source->realizationIdentity = source_realization_identity(source);
    source->realizationRoleFlags = SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL;
    source->stale = FALSE;
    source->staleReason = SoBRLDatabaseSource::STALE_NONE;
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

static int
remove_external_primary_children_except(SoBRLDatabaseSource *source,
					SoNode *keepVList,
					SoNode *keepMesh)
{
    if (!source)
	return 0;

    int removed = 0;
    for (int i = source->getNumChildren() - 1; i >= 0; i--) {
	SoNode *child = source->getChild(i);
	if (child == keepVList || child == keepMesh ||
	    node_is_source_placement_transform(child) ||
	    node_is_auxiliary_vlist(child) ||
	    node_is_auxiliary_source(child))
	    continue;
	source->removeChild(i);
	removed++;
    }
    return removed;
}

static void
clear_external_primary_vlist(SoBRLDatabaseSource *source,
			     SoBRLVListShape *shape)
{
    if (!source || !shape)
	return;

    const SbString sourceType = shape->sourceType.getValue();
    const SbString geometryKind = shape->geometryKind.getValue();
    assign_external_primary_identity(shape, source,
				     external_string_or_default(sourceType.getString(), "line-set"),
				     external_string_or_default(geometryKind.getString(), "line"));
    shape->setLineSet(NULL, NULL, 0);
    shape->setPrecisePoints(NULL, 0);
}

static void
clear_external_primary_mesh(SoBRLDatabaseSource *source,
			    SoBRLMeshShape *shape)
{
    if (!source || !shape)
	return;

    const SbString sourceType = shape->sourceType.getValue();
    const SbString geometryKind = shape->geometryKind.getValue();
    assign_external_primary_identity(shape, source,
				     external_string_or_default(sourceType.getString(), "indexed-face-set"),
				     external_string_or_default(geometryKind.getString(), "surface"));
    shape->setIndexedTriangles(NULL, 0, NULL, 0);
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
    this->discardCompactInstanceHistory();
    this->clearCompiledAssembly();
    const int before = this->getNumChildren();
    if (preserveAuxiliary) {
	remove_non_auxiliary_children(this);
    } else {
	for (int i = this->getNumChildren() - 1; i >= 0; i--) {
	    if (node_is_source_placement_transform(this->getChild(i)))
		continue;
	    this->removeChild(i);
	}
    }

    return before != this->getNumChildren() ? 1 : 0;
}

int
SoBRLDatabaseSource::clearExternalPrimaryGeometry(void)
{
    SoBRLVListShape *primaryVList = first_direct_primary_vlist_child(this);
    SoBRLMeshShape *primaryMesh = first_direct_primary_mesh_child(this);

    (void)remove_external_primary_children_except(this, primaryVList,
	    primaryMesh);

    if (!primaryVList && !primaryMesh) {
	primaryVList = new SoBRLVListShape;
	database_source_add_realized_child(this, primaryVList);
    }

    clear_external_primary_vlist(this, primaryVList);
    clear_external_primary_mesh(this, primaryMesh);
    this->clearSourceBounds();
    mark_external_primary_published_current(this);
    this->syncRealizedShapeOwnerState();
    return 1;
}

int
SoBRLDatabaseSource::publishExternalLineSet(
    const BObolExternalLineSet &lineSet)
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
	return this->clearExternalPrimaryGeometry();

    (void)this->clearRealizedGeometry(TRUE);
    SoBRLVListShape *shape = new SoBRLVListShape;
    database_source_add_realized_child(this, shape);

    assign_external_primary_identity(shape, this,
				     external_string_or_default(lineSet.sourceType, "line-set"),
				     external_string_or_default(lineSet.geometryKind, "line"));
    shape->setLineSet(lineSet.points, commands, lineSet.count);
    shape->setPrecisePoints(lineSet.precisePoints, lineSet.count);
    set_external_bounds_from_points(this, lineSet.points, lineSet.count);
    mark_external_primary_published_current(this);
    this->syncRealizedShapeOwnerState();
    return 1;
}

int
SoBRLDatabaseSource::publishExternalPointSet(
    const BObolExternalPointSet &pointSet)
{
    if (pointSet.count < 0 || (pointSet.count > 0 && !pointSet.points))
	return 0;

    if (pointSet.count == 0)
	return this->clearExternalPrimaryGeometry();

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
    return this->publishExternalLineSet(lineSet);
}

int
SoBRLDatabaseSource::publishExternalTriangleMesh(
    const BObolExternalTriangleMesh &triangleMesh)
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
	return this->clearExternalPrimaryGeometry();

    for (int i = 0; i < triangleMesh.indexCount; i++) {
	const int32_t index = triangleMesh.indices[i];
	if (index < 0 || index >= triangleMesh.pointCount)
	    return 0;
    }

    (void)this->clearRealizedGeometry(TRUE);
    SoBRLMeshShape *shape = triangleMesh.lodBacked ?
			    new SoBRLLodMeshShape : new SoBRLMeshShape;
    database_source_add_realized_child(this, shape);
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
    set_external_bounds_from_points(this, triangleMesh.points,
				    triangleMesh.pointCount);
    mark_external_primary_published_current(this);
    this->syncRealizedShapeOwnerState();
    return 1;
}

int
SoBRLDatabaseSource::publishExternalAnnotation(
    const BObolExternalAnnotation &annotation)
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

    (void)this->clearRealizedGeometry(TRUE);
    SoBRLVListShape *shape = new SoBRLVListShape;
    database_source_add_realized_child(this, shape);

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

    if (annotation.linePointCount > 0) {
	set_external_bounds_from_points(this, annotation.linePoints,
					annotation.linePointCount);
    } else {
	this->clearSourceBounds();
    }
    mark_external_primary_published_current(this);
    this->syncRealizedShapeOwnerState();
    return 1;
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

    if (realization->line_count > static_cast<size_t>(INT_MAX)) {
	primitive_realization_line_set_free(realization);
	return NULL;
    }

    if (realization->line_count == 0) {
	primitive_realization_line_set_free(realization);
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
	    primitive_realization_line_set_free(realization);
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
    shape->setLineSet(points.empty() ? NULL : points.data(),
		      commands.empty() ? NULL : commands.data(),
		      static_cast<int>(points.size()));
    shape->setPrecisePoints(precisePoints.empty() ? NULL :
			    precisePoints.data(),
			    static_cast<int>(points.size()));
    shape->geometryKind = geometryKind && geometryKind[0] ?
			  geometryKind : "line";
    primitive_realization_line_set_free(realization);
    return shape;
}

static int
publish_primitive_realization_line_set(
    SoBRLDatabaseSource *source,
    struct rt_primitive_lod_realization *realization,
    const char *sourceType)
{
    if (!source || !realization || !realization->has_line_set)
	return 0;

    if (realization->line_count == 0) {
	const int cleared = source->clearExternalPrimaryGeometry();
	primitive_realization_line_set_free(realization);
	return cleared > 0 ? 1 : 0;
    }

    SoBRLVListShape *shape =
	vlist_from_primitive_realization_line_set(realization, "line");
    if (!shape)
	return -1;
    shape->ref();

    const SoBRLVListShape *geom = shape->getGeometrySource();
    BObolExternalLineSet lineSet;
    lineSet.points = geom && geom->point.getNum() > 0 ?
		     geom->point.getValues(0) : NULL;
    lineSet.commands = geom && geom->command.getNum() > 0 ?
		       geom->command.getValues(0) : NULL;
    lineSet.precisePoints = NULL;
    lineSet.count = geom ? geom->point.getNum() : 0;
    lineSet.sourceType = sourceType && sourceType[0] ? sourceType :
			 "primitive-wireframe";
    lineSet.geometryKind = "line";
    const int published = source->publishExternalLineSet(lineSet);
    shape->unref();
    return published > 0 ? 1 : 0;
}

static void
set_external_bounds_from_vlist_shape(
    SoBRLDatabaseSource *source,
    const SoBRLVListShape *shape)
{
    if (!source || !shape)
	return;

    const SoBRLVListShape *geometrySource = shape->getGeometrySource();
    if (!geometrySource) {
	source->clearSourceBounds();
	return;
    }

    SbBox3f bounds;
    bounds.makeEmpty();
    const int pointCount = geometrySource->point.getNum();
    for (int i = 0; i < pointCount; i++)
	bounds.extendBy(geometrySource->point[i]);

    if (bounds.isEmpty()) {
	source->clearSourceBounds();
	return;
    }

    source->setSourceBoundsState(TRUE, bounds.getMin(), bounds.getMax());
}

struct primitive_submodel_publish_ctx {
    SoBRLDatabaseSource *source;
    size_t childCount;
    int failed;
};

static union tree *
primitive_submodel_wireframe_leaf(struct db_tree_state *tsp,
				  const struct db_full_path *pathp,
				  struct directory *dp,
				  void *clientData)
{
    struct primitive_submodel_publish_ctx *ctx =
	static_cast<struct primitive_submodel_publish_ctx *>(clientData);
    if (!ctx || !ctx->source || ctx->failed || !tsp || !tsp->ts_dbip ||
	!dp)
	return TREE_NULL;

    char *pathName = pathp && pathp->fp_len > 0 ?
		     db_path_to_string(pathp) : NULL;
    const char *name = (pathName && pathName[0]) ? pathName :
		       (dp->d_namep ? dp->d_namep : "submodel_leaf");

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    int haveIntern = 0;
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;

    if (rt_db_get_internal(&intern, dp, tsp->ts_dbip, NULL) < 0) {
	ctx->failed = 1;
	goto cleanup;
    }
    haveIntern = 1;

    if (!internal_payload_magic_valid(&intern) ||
	plot_internal_to_vlist_geometry(points, commands, &intern,
	    tsp->ts_ttol, tsp->ts_tol) < 0) {
	ctx->failed = 1;
	goto cleanup;
    }

    if (points.empty() || points.size() != commands.size() ||
	points.size() > static_cast<size_t>(INT_MAX)) {
	ctx->failed = 1;
	goto cleanup;
    }

    for (size_t i = 0; i < points.size(); i++) {
	point_t sourcePoint;
	point_t transformedPoint;
	VSET(sourcePoint, points[i][X], points[i][Y], points[i][Z]);
	MAT4X3PNT(transformedPoint, tsp->ts_mat, sourcePoint);
	points[i].setValue(static_cast<float>(transformedPoint[X]),
			   static_cast<float>(transformedPoint[Y]),
			   static_cast<float>(transformedPoint[Z]));
    }

    if (ctx->source->setAuxiliaryLineSet(
	    name,
	    points.empty() ? NULL : points.data(),
	    commands.empty() ? NULL : commands.data(),
	    static_cast<int>(points.size())) <= 0) {
	ctx->failed = 1;
	goto cleanup;
    }

    ctx->childCount++;

cleanup:
    if (haveIntern)
	rt_db_free_internal(&intern);
    if (pathName)
	bu_free(pathName, "BObol submodel leaf path string");
    return ctx->failed ? TREE_NULL : make_nop_tree();
}

static int
publish_primitive_submodel_wireframe(
    SoBRLDatabaseSource *source,
    struct rt_db_internal *intern,
    const struct bg_tess_tol *ttol,
    const struct bn_tol *tol)
{
    if (!source || !intern || intern->idb_type != ID_SUBMODEL ||
	!intern->idb_ptr)
	return 0;

    struct rt_submodel_internal *submodel =
	static_cast<struct rt_submodel_internal *>(intern->idb_ptr);
    RT_SUBMODEL_CK_MAGIC(submodel);

    struct db_i *dbip = DBI_NULL;
    int closeDb = 0;
    if (bu_vls_strlen(&submodel->file) != 0) {
	dbip = db_open(bu_vls_addr(&submodel->file), DB_OPEN_READONLY);
	if (dbip == DBI_NULL)
	    return -1;
	closeDb = 1;
	if (!db_is_directory_non_empty(dbip) && db_dirbuild(dbip) < 0) {
	    db_close(dbip);
	    return -1;
	}
    } else {
	RT_CK_DBI(submodel->dbip);
	dbip = const_cast<struct db_i *>(submodel->dbip);
    }

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

    (void)source->clearExternalPrimaryGeometry();
    (void)source->clearAuxiliaryShapes();

    struct db_tree_state state;
    RT_DBTS_INIT(&state);
    state.ts_dbip = dbip;
    state.ts_ttol = useTtol;
    state.ts_tol = useTol;
    MAT_COPY(state.ts_mat, submodel->root2leaf);

    struct primitive_submodel_publish_ctx ctx;
    ctx.source = source;
    ctx.childCount = 0;
    ctx.failed = 0;

    const char *argv[2];
    argv[0] = bu_vls_addr(&submodel->treetop);
    argv[1] = NULL;
    int ret = db_walk_tree_leaf_instances(dbip, 1, argv, 1, &state, 0, NULL,
	    primitive_submodel_wireframe_leaf, &ctx);

    if (closeDb)
	db_close(dbip);
    if (ret < 0 || ctx.failed)
	return -1;
    return ctx.childCount ? 1 : 0;
}

int
SoBRLDatabaseSource::publishPrimitiveWireframe(
    struct rt_db_internal *intern,
    const struct bg_tess_tol *ttol,
    const struct bn_tol *tol)
{
    if (!intern)
	return 0;
    if (!internal_payload_magic_valid(intern))
	return -1;

    if (intern->idb_type == ID_SUBMODEL)
	return publish_primitive_submodel_wireframe(this, intern, ttol, tol);

    if (intern->idb_meth && intern->idb_meth->ft_wireframe_line_set) {
	struct rt_primitive_lod_realization realization;
	memset(&realization, 0, sizeof(realization));
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
	    primitive_realization_line_set_free(&realization);
	    return -1;
	}
	return publish_primitive_realization_line_set(this, &realization,
		"line-set");
    }

    SoBRLVListShape *shape =
	vlist_from_plot_internal(intern, this, ttol, tol);
    if (!shape)
	return -1;

    (void)this->clearRealizedGeometry(TRUE);
    database_source_add_realized_child(this, shape);
    const bool annotation = primitive_is_annotation(intern->idb_type,
			    primitive_type_label(intern));
    assign_external_primary_identity(shape, this,
				     annotation ? "annotation" : "line-set",
				     annotation ? "annotation" : "line");
    set_external_bounds_from_vlist_shape(this, shape);
    mark_external_primary_published_current(this);
    this->syncRealizedShapeOwnerState();
    return 1;
}

void
SoBRLDatabaseSource::syncCompactInstanceDisplayState(void)
{
    this->rebuildCompactInstanceDisplayState(TRUE);
}

static Obol::InstanceStyle
compact_entry_style_from_source(const SoBRLDatabaseSource *source,
	const BObolCompactInstanceEntry &entry, SbBool selected,
	SbBool highlighted)
{
    Obol::InstanceStyle style;
    if (!source)
	return style;

    const BObolRealizedShapeSummary &summary = entry.shapeSummary;
    style.hasColorOverride = true;
    cad_shape_color(selected, source->selectedColor.getValue(), highlighted,
	source->highlightedColor.getValue(), summary.ghosted,
	source->ghostedColor.getValue(), source->colorOverride.getValue(),
	source->color.getValue(), entry.semantic.materialColorValid,
	entry.semantic.materialColor, source->color.getValue(),
	source->transparency.getValue(), style.color);
    style.lineWidth = source->lineWidth.getValue() > 0 ?
	static_cast<float>(source->lineWidth.getValue()) : 1.0f;
    if (source->lineStyle.getValue() != 0)
	style.linePattern = 0xcf33u;
    return style;
}

static void
compact_sync_entry_from_source(BObolCompactInstanceEntry &entry,
	const SoBRLDatabaseSource *source)
{
    if (!source)
	return;

    /* A compact entry owns occurrence metadata.  The aggregate source's
     * metadata describes the root draw request, not every leaf under it. */
    if (source->materialColorValid.getValue() &&
	(source->materialPolicy.getValue() !=
	 SoBRLDatabaseSource::MATERIAL_DATABASE ||
	 !entry.semantic.materialColorValid)) {
	entry.semantic.materialColorValid = TRUE;
	entry.semantic.materialColor = source->materialColor.getValue();
    }

    BObolRealizedShapeSummary &summary = entry.shapeSummary;
    summary.displayName = source->displayName.getValue();
    summary.drawMode = source_record_draw_mode(source);
    summary.hiddenLine = summary.drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    summary.materialRevision = source->materialRevision.getValue();
    summary.colorOverride = source->colorOverride.getValue();
    summary.color = source->color.getValue();

    entry.normalStyle = compact_entry_style_from_source(source, entry,
	FALSE, FALSE);
    entry.selectedStyle = compact_entry_style_from_source(source, entry,
	TRUE, FALSE);
    entry.highlightedStyle = compact_entry_style_from_source(source, entry,
	FALSE, TRUE);
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
	const Obol::InstanceStyle previousNormalStyle = entry.normalStyle;
	const Obol::InstanceStyle previousSelectedStyle = entry.selectedStyle;
	const Obol::InstanceStyle previousHighlightedStyle =
	    entry.highlightedStyle;
	if (syncSourceState)
	    compact_sync_entry_from_source(entry, this);
	entry.semantic.sourceInstanceKey = compact_instance_identity(entry);

	if (entry.visible != previousVisible ||
	    entry.selectable != previousSelectable)
	    entry.visibilityRevision = compact_next_revision(
		entry.visibilityRevision);
	if (entry.selected != previousSelected ||
	    entry.highlighted != previousHighlighted)
	    entry.selectionRevision = compact_next_revision(
		entry.selectionRevision);
	if (!compact_style_equal(entry.normalStyle, previousNormalStyle) ||
	    !compact_style_equal(entry.selectedStyle, previousSelectedStyle) ||
	    !compact_style_equal(entry.highlightedStyle,
		previousHighlightedStyle))
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

	if (i < this->d->compactIndex->instances.size() &&
	    this->d->compactIndex->instances[i].instance == entry.instance) {
	    this->d->compactIndex->instances[i].record.style = entry.style;
	} else {
	    auto found = std::find_if(this->d->compactIndex->instances.begin(),
		this->d->compactIndex->instances.end(),
		[&entry](const Obol::InstanceUpdate &update) {
		    return update.instance == entry.instance;
		});
	    if (found != this->d->compactIndex->instances.end())
		found->record.style = entry.style;
	}
    }
}

void
SoBRLDatabaseSource::syncCompactInstancePlacementState(void)
{
    if (!this->d->compactIndex)
	return;

    for (BObolCompactInstanceEntry &entry :
	 this->d->compactIndex->entries) {
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
	for (Obol::InstanceUpdate &update : this->d->compactIndex->instances) {
	    if (update.instance == entry.instance) {
		update.record.localToRoot = entry.localToSource;
		break;
	    }
	}
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


int
SoBRLDatabaseSource::setAuxiliaryLineSet(const char *name,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BObolAuxiliaryLineSetDisplayState *displayState)
{
    if (!name || !name[0] || count < 0 || (count > 0 && !points))
	return 0;

    SoBRLVListShape *shape = this->findAuxiliaryVListShape(name);
    if (count == 0) {
	if (!shape)
	    return 0;
	for (int i = 0; i < this->getNumChildren(); i++) {
	    if (this->getChild(i) == shape) {
		this->removeChild(i);
		return 1;
	    }
	}
	return 0;
    }

    if (!shape) {
	shape = new SoBRLVListShape;
	shape->setName(SbName(name));
	database_source_add_realized_child(this, shape);
    }

    std::vector<int32_t> fallbackCommands;
    if (!commands) {
	fallbackCommands.reserve(count);
	for (int i = 0; i < count; i++)
	    fallbackCommands.push_back(i == 0 ? SoBRLVListShape::MOVE :
				       SoBRLVListShape::DRAW);
	commands = fallbackCommands.data();
    }

    const char *sourcePath = this->path.getValue().getString();
    const uint32_t revision = this->sourceRevision.getValue();
    SbString identity = source_record_identity(this, sourcePath);
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
    shape->drawMode = source_record_draw_mode(this);
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
	sync_shape_display_state(shape, this);
    }
    sync_shape_placement_state(shape, this);
    shape->setLineSet(points, commands, count);
    sync_shape_owner_state(shape, this);
    return 1;
}


int
SoBRLDatabaseSource::setAuxiliarySourceLineSet(const char *sourcePath,
	const char *auxDisplayName,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BObolAuxiliaryLineSetDisplayState *displayState)
{
    if (!sourcePath || !sourcePath[0] || count < 0 ||
	(count > 0 && !points))
	return 0;

    SoBRLDatabaseSource *source = this->findAuxiliarySource(sourcePath);
    if (count == 0) {
	if (!source)
	    return 0;
	for (int i = 0; i < this->getNumChildren(); i++) {
	    if (this->getChild(i) == source) {
		this->removeChild(i);
		return 1;
	    }
	}
	return 0;
    }

    const std::string sourceNameStorage =
	stable_name_from_path(sourcePath, 1);
    const char *sourceName = sourceNameStorage.empty() ? sourcePath :
			     sourceNameStorage.c_str();

    if (!source) {
	source = new SoBRLDatabaseSource;
	source->setName(SbName(sourceName));
	database_source_add_realized_child(this, source);
    }

    const uint32_t revision = this->sourceRevision.getValue();
    source->configureDatabaseSourceInstance(sourcePath, sourcePath,
					    this->d->dbip, this->drawMode.getValue(), revision);
    source->auxiliarySource = TRUE;
    source->displayName = (auxDisplayName && auxDisplayName[0]) ?
			  auxDisplayName : sourceName;
    source->visible = this->visible.getValue();
    source->selected = this->selected.getValue();
    source->highlighted = this->highlighted.getValue();
    source->lineStyle = this->lineStyle.getValue();
    source->lineWidth = this->lineWidth.getValue();
    source->transparency = this->transparency.getValue();
    source->materialColorValid = this->materialColorValid.getValue();
    source->materialColor = this->materialColor.getValue();
    source->materialRevision = this->materialRevision.getValue();
    source->databaseMetadataValid = this->databaseMetadataValid.getValue();
    source->databaseRegionId = this->databaseRegionId.getValue();
    source->databaseAirCode = this->databaseAirCode.getValue();
    source->databaseMaterialId = this->databaseMaterialId.getValue();
    source->databaseLos = this->databaseLos.getValue();
    source->databaseMaterialColorValid =
	this->databaseMaterialColorValid.getValue();
    source->databaseMaterialColor = this->databaseMaterialColor.getValue();
    source->databaseMaterialShader =
	this->databaseMaterialShader.getValue();
    source->colorOverride = this->colorOverride.getValue();
    source->color = this->color.getValue();
    source->selectedColor = this->selectedColor.getValue();
    source->highlightedColor = this->highlightedColor.getValue();
    source->ghostedColor = this->ghostedColor.getValue();
    source->drawMatrixValid = this->drawMatrixValid.getValue();
    source->drawMatrix = this->drawMatrix.getValue();
    source->drawCenterValid = this->drawCenterValid.getValue();
    source->drawCenter = this->drawCenter.getValue();
    source->drawSizeValid = this->drawSizeValid.getValue();
    source->drawSize = this->drawSize.getValue();
    (void)remove_source_placement_transform(source);

    const char *shapeName = (auxDisplayName && auxDisplayName[0]) ?
			    auxDisplayName : sourceName;
    const int changed = source->setAuxiliaryLineSet(shapeName, points,
			commands, count, displayState);
    if (changed > 0) {
	source->realizedRevision = revision;
	source->realizedSourceRevision = revision;
	source->realizedInputsRevision = source->inputsRevision.getValue();
	source->realizedViewRevision = source->viewRevision.getValue();
	source->realizationStatus = REALIZED;
	source->realizationDiagnostic = "";
	source->realizationIdentity = source_realization_identity(source);
	source->stale = FALSE;
	source->staleReason = STALE_NONE;
	source->syncRealizedShapeOwnerState();
    }
    return changed;
}


int
SoBRLDatabaseSource::clearAuxiliaryShapes(void)
{
    int removed = 0;
    for (int i = this->getNumChildren() - 1; i >= 0; i--) {
	SoNode *node = this->getChild(i);
	if (!node)
	    continue;

	if (node_is_auxiliary_source(node)) {
	    this->removeChild(i);
	    removed++;
	    continue;
	}

	if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	    SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	    if (!vlist_shape_is_auxiliary(shape))
		continue;

	    this->removeChild(i);
	    removed++;
	}
    }
    return removed;
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

SbBool
SoBRLDatabaseSource::getDisplayMeshLodChangedEntries(
    uint64_t revision, std::vector<size_t> &entryIndices) const
{
    entryIndices.clear();
    if (revision == this->d->displayMeshLodRevision)
	return TRUE;
    if (!revision || revision > this->d->displayMeshLodRevision ||
	revision < this->d->displayMeshLodDeltaFloorRevision)
	return FALSE;

    for (const Impl::DisplayMeshLodDelta &delta :
	this->d->displayMeshLodDeltas) {
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

SbBool
SoBRLDatabaseSource::hasCompactInstanceIndex(void) const
{
    return (this->d->compactIndexActive && this->d->compactIndex &&
	    !this->d->compactIndex->entries.empty()) ? TRUE : FALSE;
}

SbBool
SoBRLDatabaseSource::isCompactOccurrenceRegistry(void) const
{
    return this->d->compactIndexActive && this->d->compactOccurrenceRegistry;
}

int
SoBRLDatabaseSource::getCompactInstanceCount(void) const
{
    if (!this->hasCompactInstanceIndex())
	return 0;
    return static_cast<int>(this->d->compactIndex->entries.size());
}

int
SoBRLDatabaseSource::getCompactSelectedInstanceCount(void) const
{
    if (!this->hasCompactInstanceIndex())
	return 0;
    return static_cast<int>(
	this->d->compactIndex->selectedInstances.size());
}

size_t
SoBRLDatabaseSource::getCompactExpectedInstanceCount(void) const
{
    const size_t current = this->d->compactIndex ?
	this->d->compactIndex->entries.size() : 0;
    return std::max(current, this->d->compactExpectedInstanceCount);
}

int
SoBRLDatabaseSource::getCompactPartCount(void) const
{
    if (!this->hasCompactInstanceIndex())
	return 0;
    return static_cast<int>(this->d->compactIndex->parts.size());
}

SbBool
SoBRLDatabaseSource::getCompactInstanceHandle(
    int index, BObolCompactInstanceHandle &handle) const
{
    handle = BObolCompactInstanceHandle();
    if (!this->d->compactIndex || index < 0 ||
	static_cast<size_t>(index) >= this->d->compactIndex->entries.size())
	return FALSE;

    const BObolCompactInstanceEntry &entry =
	this->d->compactIndex->entries[static_cast<size_t>(index)];
    handle.sourceNodeId = this->d->compactHandleSourceId;
    handle.instanceWord0 = entry.instance.w0;
    handle.instanceWord1 = entry.instance.w1;
    return handle.isValid();
}

SbBool
SoBRLDatabaseSource::getCompactOccurrence(
    int index, BObolCompactOccurrence &occurrence) const
{
    occurrence = BObolCompactOccurrence();
    if (!this->d->compactIndex || index < 0 ||
	static_cast<size_t>(index) >= this->d->compactIndex->entries.size())
	return FALSE;

    const BObolCompactInstanceEntry &entry =
	this->d->compactIndex->entries[static_cast<size_t>(index)];
    occurrence.geometry = entry.geometry;
    occurrence.summary = entry.shapeSummary;
    occurrence.geometryTransform = entry.geometryTransform;
    occurrence.localTransform = entry.placementTransform;
    occurrence.lodBacked = entry.lodBacked;
    occurrence.sourceMeshRequestValid = entry.sourceMeshRequestValid;
    occurrence.sourceMeshRequest = entry.sourceMeshRequest;
    occurrence.occurrenceIndex = entry.occurrenceIndex;
    occurrence.booleanOperation = entry.booleanOperation;
    return occurrence.geometry ? TRUE : FALSE;
}

SbBool
SoBRLDatabaseSource::copyCompactWireGeometry(
    std::vector<SbVec3f> &points, std::vector<int32_t> &commands) const
{
    points.clear();
    commands.clear();
    if (!this->d->compactIndex)
	return FALSE;

    for (const BObolCompactInstanceEntry &entry :
	 this->d->compactIndex->entries) {
	if (!entry.visible || !entry.geometry || !entry.geometry->wire)
	    continue;

	const Obol::WireRep &wire = *entry.geometry->wire;
	const auto appendPoint = [&entry, &points, &commands](
	    const SbVec3f &point, int32_t command) {
	    SbVec3f transformed;
	    entry.localToSource.multVecMatrix(point, transformed);
	    points.push_back(transformed);
	    commands.push_back(command);
	};

	for (size_t i = 1; i < wire.segmentPoints.size(); i += 2) {
	    appendPoint(wire.segmentPoints[i - 1], 0);
	    appendPoint(wire.segmentPoints[i], 1);
	}
	for (const Obol::WirePolyline &polyline : wire.polylines) {
	    for (size_t i = 0; i < polyline.points.size(); i++)
		appendPoint(polyline.points[i], i == 0 ? 0 : 1);
	}
    }

    return points.empty() ? FALSE : TRUE;
}

const BObolCompactInstanceEntry *
SoBRLDatabaseSource::findCompactInstanceEntry(
	const BObolCompactInstanceHandle &handle)
    const
{
    if (!this->d->compactIndex || !handle.isValid() ||
	handle.sourceNodeId != this->d->compactHandleSourceId)
	return NULL;
    Obol::InstanceId instance;
    instance.w0 = handle.instanceWord0;
    instance.w1 = handle.instanceWord1;
    const auto found = this->d->compactIndex->entryIndex.find(instance);
    if (found == this->d->compactIndex->entryIndex.end() ||
	found->second >= this->d->compactIndex->entries.size())
	return NULL;
    return &this->d->compactIndex->entries[found->second];
}

SbBool
SoBRLDatabaseSource::isCompactInstanceHandleValid(
    const BObolCompactInstanceHandle &handle) const
{
    return this->findCompactInstanceEntry(handle) ? TRUE : FALSE;
}

SbBool
SoBRLDatabaseSource::hasCompactInstanceKey(const char *occurrenceKey) const
{
    if (!this->d->compactIndex || !occurrenceKey || !occurrenceKey[0])
	return FALSE;

    const auto found =
	this->d->compactIndex->entryIndexByKey.find(occurrenceKey);
    if (found == this->d->compactIndex->entryIndexByKey.end() ||
	found->second >= this->d->compactIndex->entries.size())
	return FALSE;

    /* Keep the explicit string comparison as a defensive consistency check:
     * modern compact occurrence keys are labels derived from the instance ID,
     * not strings whose hash is the instance ID. */
    return bu_strcmp(compact_instance_identity(
	this->d->compactIndex->entries[found->second]).getString(),
	occurrenceKey) == 0 ? TRUE : FALSE;
}

SbBool
SoBRLDatabaseSource::getCompactInstanceIndex(
    const char *occurrenceKey, size_t &entryIndex) const
{
    entryIndex = 0;
    if (!this->d->compactIndex || !occurrenceKey || !occurrenceKey[0])
	return FALSE;
    const auto found =
	this->d->compactIndex->entryIndexByKey.find(occurrenceKey);
    if (found == this->d->compactIndex->entryIndexByKey.end() ||
	found->second >= this->d->compactIndex->entries.size())
	return FALSE;
    entryIndex = found->second;
    return TRUE;
}

uint64_t
SoBRLDatabaseSource::getCompactSourceRoutingId(void) const
{
    return this->d->routingId;
}

SbBool
SoBRLDatabaseSource::getCompactInstanceSummary(
    const BObolCompactInstanceHandle &handle,
    BObolCompactInstanceSummary &summary) const
{
    summary = BObolCompactInstanceSummary();
    const BObolCompactInstanceEntry *entry =
	this->findCompactInstanceEntry(handle);
    if (!entry)
	return FALSE;

    summary.valid = TRUE;
    summary.handle = handle;
    summary.path = entry->semantic.path;
    summary.sourceName = entry->semantic.sourceName;
    summary.sourceInstanceKey = compact_instance_identity(*entry);
    summary.geometryKind = entry->shapeSummary.geometryKind;
    if (entry->sourceMeshRequestValid) {
	summary.meshAssetPath = entry->sourceMeshRequest.meshAssetPath;
	summary.meshAssetName = entry->sourceMeshRequest.meshAssetName;
	summary.meshAssetBounds = entry->sourceMeshRequest.meshAssetBounds;
	summary.sourceContentHash =
	    entry->sourceMeshRequest.meshAssetContentHash;
	summary.sourceFaceCount = entry->sourceMeshRequest.faceCount;
	summary.sourcePointCount = entry->sourceMeshRequest.pointCount;
    }
    summary.localToSource = entry->localToSource;
    summary.geometryIdentity = entry->part.w0 ^
	(entry->part.w1 + 0x9e3779b97f4a7c15ULL +
	 (entry->part.w0 << 6) + (entry->part.w0 >> 2));
    summary.geometryRevision = entry->geometryRevision;
    summary.appearanceRevision = entry->appearanceRevision;
    summary.placementRevision = entry->placementRevision;
    summary.visibilityRevision = entry->visibilityRevision;
    summary.selectionRevision = entry->selectionRevision;
    summary.occurrenceIndex = entry->occurrenceIndex;
    summary.booleanOperation = entry->booleanOperation;
    summary.regionId = entry->semantic.regionId;
    summary.airCode = entry->semantic.airCode;
    summary.materialId = entry->semantic.materialId;
    summary.los = entry->semantic.los;
    summary.materialColorValid = entry->semantic.materialColorValid;
    summary.materialColor = entry->semantic.materialColor;
    summary.materialShader = entry->semantic.materialShader;
    summary.appearanceColorValid = entry->style.hasColorOverride ? TRUE : FALSE;
    summary.appearanceColor = SbColor(entry->style.color[0],
	entry->style.color[1], entry->style.color[2]);
    summary.lineStyle = entry->style.linePattern == 0xffffu ? 0 : 1;
    summary.lineWidth = entry->style.lineWidth > 0.0f ?
	static_cast<int>(entry->style.lineWidth + 0.5f) : 0;
    summary.transparency = 1.0f - entry->style.color[3];
    if (summary.transparency < 0.0f)
	summary.transparency = 0.0f;
    else if (summary.transparency > 1.0f)
	summary.transparency = 1.0f;
    summary.wireGeometry = entry->wireGeometry;
    summary.pointGeometry = entry->pointGeometry;
    summary.meshGeometry = entry->meshGeometry;
    summary.lodBacked = entry->lodBacked;
    summary.sourceMeshRequestValid = entry->sourceMeshRequestValid;
    summary.localBounds = compact_part_geometry_bounds(entry->geometry);
    summary.visible = entry->visible;
    summary.selectable = entry->selectable;
    summary.selected = entry->selected;
    summary.highlighted = entry->highlighted;
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getCompactLodInstanceSummary(
    int index, BObolCompactLodInstanceSummary &summary) const
{
    summary = BObolCompactLodInstanceSummary();
    if (!this->d->compactIndex || index < 0 ||
	static_cast<size_t>(index) >= this->d->compactIndex->entries.size())
	return FALSE;

    const BObolCompactInstanceEntry &entry =
	this->d->compactIndex->entries[static_cast<size_t>(index)];
    summary.valid = TRUE;
    summary.path = entry.semantic.path;
    summary.sourceName = entry.semantic.sourceName;
    summary.sourceInstanceKey = compact_instance_identity(entry);
    if (entry.sourceMeshRequestValid) {
	summary.meshAssetPath = entry.sourceMeshRequest.meshAssetPath;
	summary.meshAssetName = entry.sourceMeshRequest.meshAssetName;
	summary.meshAssetBounds = entry.sourceMeshRequest.meshAssetBounds;
	summary.sourceContentHash =
	    entry.sourceMeshRequest.meshAssetContentHash;
	summary.sourceFaceCount = entry.sourceMeshRequest.faceCount;
	summary.sourcePointCount = entry.sourceMeshRequest.pointCount;
    }
    summary.localToSource = entry.sourceMeshRequestValid ?
	compact_mesh_asset_matrix(this, entry) : entry.localToSource;
    summary.localBounds = entry.sourceMeshRequestValid &&
	!entry.sourceMeshRequest.meshAssetBounds.isEmpty() ?
	entry.sourceMeshRequest.meshAssetBounds :
	compact_part_geometry_bounds(entry.geometry);
    summary.meshGeometry = entry.meshGeometry;
    summary.lodBacked = entry.lodBacked;
    summary.sourceMeshRequestValid = entry.sourceMeshRequestValid;
    summary.visible = entry.visible;
    summary.selected = entry.selected;
    summary.highlighted = entry.highlighted;
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getCompactSourceMeshRequest(
    int index, BObolSourceMeshRequest &request) const
{
    request.clear();
    if (!this->d->compactIndex || index < 0 ||
	static_cast<size_t>(index) >= this->d->compactIndex->entries.size())
	return FALSE;
    const BObolCompactInstanceEntry &entry =
	this->d->compactIndex->entries[static_cast<size_t>(index)];
    if (!entry.sourceMeshRequestValid)
	return FALSE;
    request = entry.sourceMeshRequest;
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getCompactLodPlanningSummary(
    int index, BObolCompactLodPlanningSummary &summary) const
{
    summary = BObolCompactLodPlanningSummary();
    if (!this->d->compactIndex || index < 0 ||
	static_cast<size_t>(index) >= this->d->compactIndex->entries.size())
	return FALSE;

    const BObolCompactInstanceEntry &entry =
	this->d->compactIndex->entries[static_cast<size_t>(index)];
    summary.valid = TRUE;
    summary.sourceInstanceKey = compact_instance_identity(entry);
    summary.localToSource = entry.sourceMeshRequestValid ?
	compact_mesh_asset_matrix(this, entry) : entry.localToSource;
    summary.localBounds = entry.sourceMeshRequestValid &&
	!entry.sourceMeshRequest.meshAssetBounds.isEmpty() ?
	entry.sourceMeshRequest.meshAssetBounds :
	compact_part_geometry_bounds(entry.geometry);
    summary.meshGeometry = entry.meshGeometry;
    summary.lodBacked = entry.lodBacked;
    summary.sourceMeshRequestValid = entry.sourceMeshRequestValid;
    summary.visible = entry.visible;
    summary.selected = entry.selected;
    summary.highlighted = entry.highlighted;
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getCompactLodPlanningSummaryForKey(
    const char *occurrenceKey, BObolCompactLodPlanningSummary &summary) const
{
    summary = BObolCompactLodPlanningSummary();
    if (!this->d->compactIndex || !occurrenceKey || !occurrenceKey[0])
	return FALSE;
    const auto found =
	this->d->compactIndex->entryIndexByKey.find(occurrenceKey);
    if (found == this->d->compactIndex->entryIndexByKey.end() ||
	found->second >= this->d->compactIndex->entries.size())
	return FALSE;
    return this->getCompactLodPlanningSummary(
	static_cast<int>(found->second), summary);
}

SbBool
SoBRLDatabaseSource::copyCompactInstanceEditGeometry(
    const BObolCompactInstanceHandle &handle,
    std::vector<SbVec3f> &points,
    std::vector<int32_t> &commands,
    BObolCompactInstanceSummary &summary) const
{
    points.clear();
    commands.clear();
    summary = BObolCompactInstanceSummary();

    const BObolCompactInstanceEntry *entry =
	this->findCompactInstanceEntry(handle);
    if (!entry || !entry->geometry ||
	!this->getCompactInstanceSummary(handle, summary))
	return FALSE;

    const auto appendPoint = [&entry, &points, &commands](
	const SbVec3f &point, int32_t command) {
	SbVec3f transformed;
	entry->localToSource.multVecMatrix(point, transformed);
	points.push_back(transformed);
	commands.push_back(command);
    };

    if (entry->geometry->wire) {
	const Obol::WireRep &wire = *entry->geometry->wire;
	for (size_t i = 1; i < wire.segmentPoints.size(); i += 2) {
	    appendPoint(wire.segmentPoints[i - 1], 0);
	    appendPoint(wire.segmentPoints[i], 1);
	}
	for (const Obol::WirePolyline &polyline : wire.polylines) {
	    for (size_t i = 0; i < polyline.points.size(); i++)
		appendPoint(polyline.points[i], i == 0 ? 0 : 1);
	}
    }

    if (entry->geometry->points) {
	const Obol::PointRep &pointRep = *entry->geometry->points;
	for (const SbVec3f &point : pointRep.positions)
	    appendPoint(point, 2);
    }

    /* Mesh-only compact occurrences still need an editable visual.  Build a
     * transient triangle-edge preview without adding a persistent mesh shape
     * to the compact index. */
    if (points.empty() && entry->geometry->shaded) {
	const Obol::TriMesh &mesh = *entry->geometry->shaded;
	for (size_t i = 0; i + 2 < mesh.indices.size(); i += 3) {
	    const uint32_t a = mesh.indices[i];
	    const uint32_t b = mesh.indices[i + 1];
	    const uint32_t c = mesh.indices[i + 2];
	    if (a >= mesh.positions.size() || b >= mesh.positions.size() ||
		c >= mesh.positions.size())
		continue;
	    appendPoint(mesh.positions[a], 0);
	    appendPoint(mesh.positions[b], 1);
	    appendPoint(mesh.positions[c], 1);
	    appendPoint(mesh.positions[a], 1);
	}
    }

    if (points.empty()) {
	summary = BObolCompactInstanceSummary();
	return FALSE;
    }
    return TRUE;
}

template <typename Visitor>
static void
compact_visit_entries_for_path(const BObolCompactInstanceIndex *index,
	const char *queryPath, SbBool includeDescendants, Visitor visitor)
{
    if (!index)
	return;

    const char *query = database_source_skip_leading_slash(
	queryPath ? queryPath : "");
    if (!query[0]) {
	for (size_t entryIndex = 0; entryIndex < index->entries.size();
		entryIndex++)
	    visitor(entryIndex);
	return;
    }

    const bool leafQuery = !strchr(query, '/') && !strchr(query, '@');
    if (leafQuery) {
	auto leafEntries = index->entryIndicesByLeaf.find(query);
	if (leafEntries != index->entryIndicesByLeaf.end()) {
	    for (size_t entryIndex : leafEntries->second)
		visitor(entryIndex);
	}
	if (!includeDescendants)
	    return;
    }

    const std::string pathKey(query);
    const size_t prefixLength = pathKey.size();
    auto entryIt = index->entryIndexByOrderedPath.lower_bound(pathKey);
    for (; entryIt != index->entryIndexByOrderedPath.end(); ++entryIt) {
	const char *candidate = entryIt->first.c_str();
	if (bu_strncmp(candidate, pathKey.c_str(), prefixLength) != 0)
	    break;
	const char suffix = candidate[prefixLength];
	if (!includeDescendants && suffix != '\0')
	    break;
	if (includeDescendants && suffix != '\0' && suffix != '/' &&
	    suffix != '@')
	    continue;
	const size_t entryIndex = entryIt->second;
	if (entryIndex >= index->entries.size())
	    continue;
	if (leafQuery && database_source_leaf_component(
		index->entries[entryIndex].semantic.path) == pathKey)
	    continue;
	visitor(entryIndex);
	if (!includeDescendants)
	    continue;
    }
}

int
SoBRLDatabaseSource::getCompactInstanceCountForPath(const char *queryPath,
    SbBool includeDescendants) const
{
    if (!this->d->compactIndex)
	return 0;
    int count = 0;
    compact_visit_entries_for_path(this->d->compactIndex, queryPath,
	includeDescendants, [&count](size_t UNUSED(entryIndex)) {
	    count++;
	});
    return count;
}

SbBool
SoBRLDatabaseSource::getCompactInstanceBoundsForPath(const char *queryPath,
    SbBool includeDescendants, SbBox3f &bounds) const
{
    bounds.makeEmpty();
    if (!this->d->compactIndex || !this->visible.getValue())
	return FALSE;

    compact_visit_entries_for_path(this->d->compactIndex, queryPath,
	includeDescendants, [this, &bounds](size_t entryIndex) {
	const BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[entryIndex];
	if (!entry.visible)
	    return;
	const SbBox3f localBounds = compact_part_geometry_bounds(entry.geometry);
	if (!localBounds.isEmpty())
	    bounds.extendBy(database_source_transform_bounds(localBounds,
		entry.localToSource));
    });
    return bounds.isEmpty() ? FALSE : TRUE;
}

int
SoBRLDatabaseSource::setCompactInstanceDisplayStateForPath(const char *queryPath,
    SbBool includeDescendants,
    int visibleValid, SbBool nextVisible,
    int selectedValid, SbBool nextSelected,
    int highlightedValid, SbBool nextHighlighted)
{
    if (!this->d->compactIndex)
	return 0;
    int changed = 0;
    bool anyVisibilityChanged = false;
    std::vector<size_t> visibilityChangedEntries;

    compact_visit_entries_for_path(this->d->compactIndex, queryPath,
	includeDescendants, [this, visibleValid, nextVisible, selectedValid,
	nextSelected, highlightedValid, nextHighlighted, &changed,
	&anyVisibilityChanged, &visibilityChangedEntries](size_t entryIndex) {
	BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[entryIndex];
	bool visibilityChanged = false;
	bool selectionChanged = false;
	/*
	 * Retirement of the synthetic whole-target extent is monotonic within
	 * one compact source epoch.  A redraw synchronizes the root path as
	 * visible after the authoritative stream has completed; treating that
	 * user-facing root visibility as authored visibility for every internal
	 * record used to resurrect the retired overview and leave one box over
	 * the finished model.  A genuinely new realization installs a new index
	 * (and therefore a newly authored overview), so there is no valid reason
	 * to revive this internal record in place.
	 */
	const bool retiredOverview =
	    BU_STR_EQUAL(entry.shapeSummary.recordRole.getString(),
		"lod-overview") && !entry.authoredVisible;
	if (visibleValid && !(nextVisible && retiredOverview) &&
	    entry.authoredVisible != nextVisible) {
	    entry.authoredVisible = nextVisible;
	    entry.visible = nextVisible;
	    visibilityChanged = true;
	    anyVisibilityChanged = true;
	    visibilityChangedEntries.push_back(entryIndex);
	    changed++;
	}
	if (selectedValid && entry.selected != nextSelected) {
	    entry.selected = nextSelected;
	    selectionChanged = true;
	    changed++;
	}
	if (highlightedValid && entry.highlighted != nextHighlighted) {
	    entry.highlighted = nextHighlighted;
	    selectionChanged = true;
	    changed++;
	}
	if (visibilityChanged)
	    entry.visibilityRevision = compact_next_revision(
		entry.visibilityRevision);
	if (selectionChanged) {
	    entry.selectionRevision = compact_next_revision(
		entry.selectionRevision);
	    entry.style = compact_effective_style(entry);
	}
    });
    if (changed) {
	if (visibleValid && this->d->compactVisibilityFrontierActive)
	    (void)this->reapplyCompactInstanceVisibilityFrontier();
	else
	    this->rebuildCompactInstanceDisplayState(FALSE);
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
	if (anyVisibilityChanged &&
	    !this->d->compactVisibilityFrontierActive)
	    this->markDisplayMeshLodDirty(visibilityChangedEntries);
	this->touch();
    }
    return changed;
}


int
SoBRLDatabaseSource::reapplyCompactInstanceVisibilityFrontier(
    size_t firstEntry)
{
    if (!this->d->compactIndex)
	return 0;

    BObolCompactInstanceIndex &index = *this->d->compactIndex;
    firstEntry = std::min(firstEntry, index.entries.size());
    const size_t candidateCount = index.entries.size() - firstEntry;
    std::vector<SbBool> allowed(candidateCount,
	this->d->compactVisibilityFrontierActive ?
	    this->d->compactVisibilityFrontierDefault : TRUE);
    if (this->d->compactVisibilityFrontierActive) {
	for (size_t ruleIndex = 0;
	     ruleIndex < this->d->compactVisibilityFrontier.size();
	     ruleIndex++) {
	    const SbString &frontierPath =
		this->d->compactVisibilityFrontier[ruleIndex];
	    const SbBool ruleVisible =
		ruleIndex < this->d->compactVisibilityFrontierStates.size() ?
		this->d->compactVisibilityFrontierStates[ruleIndex] : TRUE;
	    compact_visit_entries_for_path(&index,
		frontierPath.getString(), TRUE,
		[&allowed, ruleVisible, firstEntry](size_t entryIndex) {
		    if (entryIndex >= firstEntry &&
			entryIndex - firstEntry < allowed.size())
			allowed[entryIndex - firstEntry] = ruleVisible;
		});
	}
    }

    int changed = 0;
    std::vector<size_t> changedEntries;
    for (size_t i = firstEntry; i < index.entries.size(); i++) {
	BObolCompactInstanceEntry &entry = index.entries[i];
	const SbBool nextVisible =
	    entry.authoredVisible && allowed[i - firstEntry];
	if (entry.visible == nextVisible)
	    continue;
	entry.visible = nextVisible;
	entry.visibilityRevision = compact_next_revision(
	    entry.visibilityRevision);
	compact_sync_shape_summary_state(entry);
	changedEntries.push_back(i);
	changed++;
    }
    if (changed) {
	/* Visibility is an instance-set delta, not a geometry or style change.
	 * Rebuilding all display state here needlessly revisits every style and
	 * instance record and can make a one-occurrence erase look like a scene
	 * rebuild to the renderer. */
	if (firstEntry == 0) {
	    index.hiddenInstances.clear();
	    for (const BObolCompactInstanceEntry &entry : index.entries) {
		if (!entry.visible)
		    index.hiddenInstances.push_back(entry.instance);
	    }
	} else {
	    /* Streamed candidates are newly appended and therefore cannot
	     * already occur in hiddenInstances.  Publish only the hidden
	     * additions rather than rescanning every earlier occurrence. */
	    for (const size_t entryIndex : changedEntries) {
		const BObolCompactInstanceEntry &entry =
		    index.entries[entryIndex];
		if (!entry.visible)
		    index.hiddenInstances.push_back(entry.instance);
	    }
	}
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
	this->markDisplayMeshLodDirty(changedEntries);
	this->touch();
    }
    return changed;
}


int
SoBRLDatabaseSource::setCompactInstanceVisibilityFrontier(
    const std::vector<SbString> &paths)
{
    bool same = this->d->compactVisibilityFrontierActive &&
	this->d->compactVisibilityFrontierDefault == FALSE &&
	this->d->compactVisibilityFrontier.size() == paths.size() &&
	this->d->compactVisibilityFrontierStates.size() == paths.size();
    for (size_t i = 0; same && i < paths.size(); i++)
	same = database_source_string_equal(
	    this->d->compactVisibilityFrontier[i], paths[i].getString()) &&
	    this->d->compactVisibilityFrontierStates[i] == TRUE;
    if (same)
	return 0;

    this->d->compactVisibilityFrontier = paths;
    this->d->compactVisibilityFrontierStates.assign(paths.size(), TRUE);
    this->d->compactVisibilityFrontierDefault = FALSE;
    this->d->compactVisibilityFrontierActive = TRUE;
    const int changed = this->reapplyCompactInstanceVisibilityFrontier();
    return changed > 0 ? changed : 1;
}


int
SoBRLDatabaseSource::setCompactInstanceVisibilityOverrides(
    const std::vector<SbString> &paths,
    const std::vector<SbBool> &states)
{
    if (paths.size() != states.size())
	return 0;

    bool same = this->d->compactVisibilityFrontierActive &&
	this->d->compactVisibilityFrontierDefault == TRUE &&
	this->d->compactVisibilityFrontier.size() == paths.size() &&
	this->d->compactVisibilityFrontierStates.size() == states.size();
    for (size_t i = 0; same && i < paths.size(); i++) {
	same = database_source_string_equal(
	    this->d->compactVisibilityFrontier[i], paths[i].getString()) &&
	    this->d->compactVisibilityFrontierStates[i] == states[i];
    }
    if (same)
	return 0;

    this->d->compactVisibilityFrontier = paths;
    this->d->compactVisibilityFrontierStates = states;
    this->d->compactVisibilityFrontierDefault = TRUE;
    this->d->compactVisibilityFrontierActive = TRUE;
    const int changed = this->reapplyCompactInstanceVisibilityFrontier();
    return changed > 0 ? changed : 1;
}


int
SoBRLDatabaseSource::clearCompactInstanceVisibilityFrontier(void)
{
    if (!this->d->compactVisibilityFrontierActive)
	return 0;
    this->d->compactVisibilityFrontierActive = FALSE;
    this->d->compactVisibilityFrontierDefault = FALSE;
    this->d->compactVisibilityFrontier.clear();
    this->d->compactVisibilityFrontierStates.clear();
    const int changed = this->reapplyCompactInstanceVisibilityFrontier();
    return changed > 0 ? changed : 1;
}


SbBool
SoBRLDatabaseSource::hasCompactInstanceVisibilityFrontier(void) const
{
    return this->d->compactVisibilityFrontierActive;
}


int
SoBRLDatabaseSource::reapplyCompactInstanceSelectedPaths(
    size_t firstEntry)
{
    if (!this->d->compactIndex)
	return 0;

    BObolCompactInstanceIndex &index = *this->d->compactIndex;
    firstEntry = std::min(firstEntry, index.entries.size());
    const size_t candidateCount = index.entries.size() - firstEntry;
    std::vector<unsigned char> touched(candidateCount, 0);
    std::vector<unsigned char> selectedMask(candidateCount, 0);

    /* Existing selected entries are the only candidates that can become
     * deselected.  New selection paths contribute candidates through the
     * sorted path index, so selecting one row does not scan every occurrence
     * in a large assembly. */
    if (firstEntry == 0) {
	for (const Obol::InstanceId &instance : index.selectedInstances) {
	    auto found = index.entryIndex.find(instance);
	    if (found != index.entryIndex.end() &&
		found->second < touched.size())
		touched[found->second] = 1;
	}
    }

    const auto selectPath =
	[&index, &touched, &selectedMask, firstEntry](
	const char *queryPath) {
	if (!queryPath || !queryPath[0])
	    return;
	compact_visit_entries_for_path(&index, queryPath, TRUE,
	    [&touched, &selectedMask, firstEntry](size_t entryIndex) {
		if (entryIndex < firstEntry ||
		    entryIndex - firstEntry >= touched.size())
		    return;
		touched[entryIndex - firstEntry] = 1;
		selectedMask[entryIndex - firstEntry] = 1;
	    });
    };
    for (const SbString &selectedPath : this->d->compactSelectedPaths) {
	const char *selectedPathString = database_source_skip_leading_slash(
	    selectedPath.getString());
	if (!selectedPathString || !selectedPathString[0])
	    continue;
	selectPath(selectedPathString);
	const std::string semantic =
	    database_source_db_path_without_instance_suffixes(
		selectedPathString);
	if (!semantic.empty() && semantic != selectedPathString)
	    selectPath(semantic.c_str());
    }
    if (getenv("BOBOL_SELECTION_DEBUG")) {
	size_t touchedCount = 0;
	size_t selectedCount = 0;
	for (size_t i = 0; i < touched.size(); i++) {
	    touchedCount += touched[i] ? 1u : 0u;
	    selectedCount += selectedMask[i] ? 1u : 0u;
	}
	bu_log("[obol-selection] source=%s entries=%zu paths=%zu "
	    "touched=%zu selected=%zu\n",
	    this->path.getValue().getString(), index.entries.size(),
	    this->d->compactSelectedPaths.size(), touchedCount,
	    selectedCount);
	for (const SbString &selectedPath : this->d->compactSelectedPaths)
	    bu_log("[obol-selection]   path=%s\n",
		selectedPath.getString());
    }

    int changed = 0;
    std::vector<size_t> changedEntries;
    for (size_t i = firstEntry; i < index.entries.size(); i++) {
	if (!touched[i - firstEntry])
	    continue;
	BObolCompactInstanceEntry &entry = index.entries[i];
	const SbBool nextSelected =
	    selectedMask[i - firstEntry] ? TRUE : FALSE;
	if (entry.selected == nextSelected)
	    continue;
	entry.selected = nextSelected;
	entry.selectionRevision = compact_next_revision(
	    entry.selectionRevision);
	entry.style = compact_effective_style(entry);
	compact_sync_shape_summary_state(entry);
	if (i < index.instances.size() &&
	    index.instances[i].instance == entry.instance) {
	    index.instances[i].record.style = entry.style;
	} else {
	    auto found = std::find_if(index.instances.begin(),
		index.instances.end(),
		[&entry](const Obol::InstanceUpdate &update) {
		    return update.instance == entry.instance;
		});
	    if (found != index.instances.end())
		found->record.style = entry.style;
	}
	changedEntries.push_back(i);
	changed++;
    }
    if (changed) {
	std::unordered_set<Obol::InstanceId, std::hash<Obol::InstanceId>>
	    changedInstances;
	changedInstances.reserve(changedEntries.size());
	for (size_t entryIndex : changedEntries)
	    changedInstances.insert(index.entries[entryIndex].instance);

	std::vector<Obol::InstanceId> nextSelected;
	nextSelected.reserve(index.selectedInstances.size() +
	    changedEntries.size());
	for (const Obol::InstanceId &instance : index.selectedInstances) {
	    if (changedInstances.find(instance) == changedInstances.end())
		nextSelected.push_back(instance);
	}
	for (size_t entryIndex : changedEntries) {
	    const BObolCompactInstanceEntry &entry = index.entries[entryIndex];
	    if (entry.selected)
		nextSelected.push_back(entry.instance);
	}
	index.selectedInstances.swap(nextSelected);
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty(changedEntries);
	this->touch();
    }
    return changed;
}


int
SoBRLDatabaseSource::syncCompactInstanceSelectedPaths(
    const std::vector<SbString> &paths)
{
    bool same = this->d->compactSelectedPaths.size() == paths.size();
    for (size_t i = 0; same && i < paths.size(); i++)
	same = database_source_string_equal(
	    this->d->compactSelectedPaths[i], paths[i].getString());
    if (same)
	return 0;

    this->d->compactSelectedPaths = paths;
    const int changed = this->reapplyCompactInstanceSelectedPaths();
    return changed > 0 ? changed : 1;
}


int
SoBRLDatabaseSource::setCompactInstanceSelectableForPath(
    const char *queryPath, SbBool includeDescendants, SbBool nextSelectable)
{
    if (!this->d->compactIndex)
	return 0;
    int changed = 0;

    compact_visit_entries_for_path(this->d->compactIndex, queryPath,
	includeDescendants, [this, nextSelectable, &changed](size_t entryIndex) {
	BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[entryIndex];
	if (entry.selectable == nextSelectable)
	    return;
	entry.selectable = nextSelectable;
	entry.selectionRevision = compact_next_revision(
	    entry.selectionRevision);
	changed++;
    });
    if (changed) {
	this->rebuildCompactInstanceDisplayState(FALSE);
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
	this->touch();
    }
    return changed;
}

int
SoBRLDatabaseSource::setCompactInstanceRegionIdForPath(const char *queryPath,
    SbBool includeDescendants, int regionId)
{
    if (!this->d->compactIndex)
	return 0;
    int changed = 0;

    compact_visit_entries_for_path(this->d->compactIndex, queryPath,
	includeDescendants, [this, regionId, &changed](size_t entryIndex) {
	BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[entryIndex];
	if (entry.semantic.regionId == regionId)
	    return;
	entry.semantic.regionId = regionId;
	entry.appearanceRevision = compact_next_revision(
	    entry.appearanceRevision);
	compact_sync_shape_summary_state(entry);
	changed++;
    });
    if (changed) {
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
    }
    return changed;
}

int
SoBRLDatabaseSource::setCompactInstanceRegionMetadataForPath(
    const char *queryPath, SbBool includeDescendants, int regionId,
    int airCode, int materialId, int los)
{
    if (!this->d->compactIndex)
	return 0;
    int changed = 0;

    compact_visit_entries_for_path(this->d->compactIndex, queryPath,
	includeDescendants, [this, regionId, airCode, materialId, los,
	&changed](size_t entryIndex) {
	BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[entryIndex];
	if (entry.semantic.regionId == regionId &&
	     entry.semantic.airCode == airCode &&
	     entry.semantic.materialId == materialId &&
	     entry.semantic.los == los)
	    return;
	entry.semantic.regionId = regionId;
	entry.semantic.airCode = airCode;
	entry.semantic.materialId = materialId;
	entry.semantic.los = los;
	entry.appearanceRevision = compact_next_revision(
	    entry.appearanceRevision);
	compact_sync_shape_summary_state(entry);
	changed++;
    });
    if (changed) {
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
    }
    return changed;
}

int
SoBRLDatabaseSource::setCompactInstanceMetadataForPath(const char *queryPath,
    SbBool includeDescendants, int regionId, int airCode, int materialId,
    int los, SbBool nextMaterialColorValid, const SbColor &nextMaterialColor,
    const SbString &materialShader)
{
    if (!this->d->compactIndex)
	return 0;
    int changed = 0;

    compact_visit_entries_for_path(this->d->compactIndex, queryPath,
	includeDescendants, [this, regionId, airCode, materialId, los,
	nextMaterialColorValid, &nextMaterialColor, &materialShader,
	&changed](size_t entryIndex) {
	BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[entryIndex];
	const bool same = entry.semantic.regionId == regionId &&
	    entry.semantic.airCode == airCode &&
	    entry.semantic.materialId == materialId &&
	    entry.semantic.los == los &&
	    entry.semantic.materialColorValid == nextMaterialColorValid &&
	    (!nextMaterialColorValid || database_source_color_equal(
		entry.semantic.materialColor, nextMaterialColor)) &&
	    bu_strcmp(entry.semantic.materialShader.getString(),
		materialShader.getString()) == 0;
	if (same)
	    return;
	entry.semantic.regionId = regionId;
	entry.semantic.airCode = airCode;
	entry.semantic.materialId = materialId;
	entry.semantic.los = los;
	entry.semantic.materialColorValid = nextMaterialColorValid;
	entry.semantic.materialColor = nextMaterialColorValid ? nextMaterialColor :
	    SbColor(1.0f, 1.0f, 1.0f);
	entry.semantic.materialShader = materialShader;
	entry.normalStyle = compact_entry_style_from_source(this, entry,
	    FALSE, FALSE);
	entry.selectedStyle = compact_entry_style_from_source(this, entry,
	    TRUE, FALSE);
	entry.highlightedStyle = compact_entry_style_from_source(this, entry,
	    FALSE, TRUE);
	entry.style = compact_effective_style(entry);
	entry.appearanceRevision = compact_next_revision(
	    entry.appearanceRevision);
	compact_sync_shape_summary_state(entry);
	changed++;
    });
    if (changed) {
	this->rebuildCompactInstanceDisplayState(FALSE);
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
    }
    return changed;
}

int
SoBRLDatabaseSource::setCompactSubtractLineStyle(int nextLineStyle)
{
    if (!this->d->compactIndex)
	return 0;
    int changed = 0;
    const uint16_t pattern = nextLineStyle != 0 ? 0xcf33u : 0xffffu;
    for (BObolCompactInstanceEntry &entry : this->d->compactIndex->entries) {
	if (entry.booleanOperation != BOOLEAN_SUBTRACT || !entry.wireGeometry ||
	    entry.normalStyle.linePattern == pattern)
	    continue;
	entry.normalStyle.linePattern = pattern;
	entry.selectedStyle.linePattern = pattern;
	entry.highlightedStyle.linePattern = pattern;
	entry.style = compact_effective_style(entry);
	entry.appearanceRevision = compact_next_revision(
	    entry.appearanceRevision);
	changed++;
    }
    if (changed) {
	this->rebuildCompactInstanceDisplayState(FALSE);
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
    }
    return changed;
}

int
SoBRLDatabaseSource::refreshCompactObjectGeometry(
    const char *objectPath, uint32_t nextSourceRevision)
{
    if (!this->isCompactOccurrenceRegistry() || !this->d->compactIndex ||
	!this->d->dbip || !objectPath || !objectPath[0])
	return 0;

    const char *objectName = strrchr(objectPath, '/');
    objectName = objectName && objectName[1] ? objectName + 1 : objectPath;
    if (!objectName[0])
	return 0;

    struct directory *dp = db_lookup(this->d->dbip, objectName, LOOKUP_QUIET);
    if (!dp)
	return -1;
    if (dp->d_flags & RT_DIR_COMB) {
	BObolDatabaseSourceRealizationCache cache;
	cache.preserveCompactSourceOnFailure = true;
	this->seedCompactRealizationCache(&cache);

	const uint32_t previousRevision = this->sourceRevision.getValue();
	const uint32_t revision = nextSourceRevision ? nextSourceRevision :
	    previousRevision + 1;
	this->detachFieldSensors();
	this->sourceRevision = revision;
	this->attachFieldSensors();

	const int realized = this->usesMeshRealization() ?
	    bobol_database_source_realize_mesh_compact_with_cache(this,
		&cache) :
	    bobol_database_source_realize_wireframe_compact_with_cache(this,
		&cache);
	if (realized <= 0) {
	    this->detachFieldSensors();
	    this->sourceRevision = previousRevision;
	    this->attachFieldSensors();
	    mark_source_realized_current(this);
	    return -1;
	}
	this->markCadBatchDirty();
	return this->getCompactInstanceCount();
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
    Obol::PartGeometry generated;
    bool geometryValid = false;
    bool directWire = false;
    bool directMesh = false;
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
		&validInternal.local, this, localBounds, generated) != 0;
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

    if (sharedWire)
	sharedWire->ref();
    if (sharedMesh)
	sharedMesh->ref();
    const uint32_t revision = nextSourceRevision ? nextSourceRevision :
	this->sourceRevision.getValue() + 1;
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
	    if (sharedWire)
		sharedWire->unref();
	    if (sharedMesh)
		sharedMesh->unref();
	    return -1;
	}
    }
    if (!geometryValid && sharedWire) {
	SoBRLVListShape *geometryShape = new SoBRLVListShape;
	geometryShape->ref();
	geometryShape->setSharedGeometry(sharedWire);
	geometryValid = cad_vlist_part_geometry(geometryShape, generated) != 0;
	geometryShape->unref();
    } else if (!geometryValid && sharedMesh) {
	SoBRLMeshShape *geometryShape = new SoBRLMeshShape;
	geometryShape->ref();
	geometryShape->setSharedGeometry(sharedMesh);
	geometryShape->hiddenLine =
	    sourceDrawMode == BOBOL_LOD_DRAW_HIDDEN_LINE ? TRUE : FALSE;
	geometryShape->drawMode = sourceDrawMode;
	geometryValid = cad_mesh_part_geometry(geometryShape, generated) != 0;
	geometryShape->unref();
    }
    if (!geometryValid) {
	if (sharedWire)
	    sharedWire->unref();
	if (sharedMesh)
	    sharedMesh->unref();
	return -1;
    }

    std::string partKey;
    const char *partKind = (sharedWire || directWire) ? "wire" : "mesh";
    if (!cad_part_key_for_geometry(partKind,
	generated, partKey)) {
	if (sharedWire)
	    sharedWire->unref();
	if (sharedMesh)
	    sharedMesh->unref();
	return -1;
    }
	Obol::PartId partId;
	std::shared_ptr<const Obol::PartGeometry> geometry;
	const std::map<std::string, Obol::PartId>::const_iterator existingPart =
	this->d->compactIndex->partIdByKey.find(partKey);
	if (existingPart != this->d->compactIndex->partIdByKey.end()) {
	    partId = existingPart->second;
	    for (const BObolCompactPartReference &partRef :
		 this->d->compactIndex->parts) {
		if (partRef.part == partId) {
		    geometry = partRef.geometry;
		    break;
		}
	    }
	}
	if (!geometry) {
	    partId = Obol::CadIdBuilder::hash128(partKey);
	    geometry = std::make_shared<const Obol::PartGeometry>(
		std::move(generated));
	    BObolCompactPartReference partRef;
	    partRef.part = partId;
	    partRef.geometry = geometry;
	    this->d->compactIndex->parts.push_back(partRef);
	}
	this->d->compactIndex->partIdByKey[partKey] = partId;
	this->d->compactIndex->partIdByGeometry[geometry.get()] = partId;

    std::unordered_set<Obol::PartId, std::hash<Obol::PartId>> oldParts;
    for (size_t index : matching) {
	BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[index];
	oldParts.insert(entry.part);
	compact_index_bounds_remove(*this->d->compactIndex, entry);
	if (entry.meshGeometry &&
	    this->d->compactIndex->shadedCount > 0)
	    this->d->compactIndex->shadedCount--;
	if ((entry.wireGeometry || entry.pointGeometry) &&
	    this->d->compactIndex->wireCount > 0)
	    this->d->compactIndex->wireCount--;
	std::unordered_map<Obol::PartId, size_t, std::hash<Obol::PartId>>::iterator
	    oldPartCount = this->d->compactIndex->partReferenceCounts.find(entry.part);
	if (oldPartCount != this->d->compactIndex->partReferenceCounts.end() &&
	    oldPartCount->second > 0)
	    oldPartCount->second--;
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
	entry.lodBacked = replacementLodBacked ? TRUE : FALSE;
	entry.shapeSummary.geometryKind =
	    replacementLodBacked ? "aabb" :
	    ((sharedWire || directWire) ? "line" : "surface");
	const SbBool oldSourceMeshRequestValid = entry.sourceMeshRequestValid;
	entry.sourceMeshRequestValid = replacementSourceMeshRequestValid ?
	    TRUE : FALSE;
	if (entry.sourceMeshRequestValid)
	    entry.sourceMeshRequest = replacementSourceMeshRequest;
	else
	    entry.sourceMeshRequest.clear();
	if (oldSourceMeshRequestValid != entry.sourceMeshRequestValid) {
	    if (entry.sourceMeshRequestValid)
		this->d->compactIndex->sourceMeshRequestCount++;
	    else if (this->d->compactIndex->sourceMeshRequestCount > 0)
		this->d->compactIndex->sourceMeshRequestCount--;
	}
	entry.geometryRevision = compact_next_revision(entry.geometryRevision);
	entry.semantic.sourceId = revision;
	compact_sync_shape_summary(entry);
	if (entry.sourceMeshRequestValid) {
	    compact_source_mesh_request_sync(entry.sourceMeshRequest,
		entry.shapeSummary);
	    compact_summary_lod_from_source_mesh_request(entry.shapeSummary,
		entry.sourceMeshRequest);
	}
	compact_index_bounds_add(*this->d->compactIndex, entry);
	if (entry.meshGeometry)
	    this->d->compactIndex->shadedCount++;
	if (entry.wireGeometry || entry.pointGeometry)
	    this->d->compactIndex->wireCount++;
	if (index < this->d->compactIndex->instances.size()) {
	    this->d->compactIndex->instances[index].record.part = partId;
	    this->d->compactIndex->instances[index].record.localToRoot =
		entry.localToSource;
	}
    }
    this->d->compactIndex->partReferenceCounts[partId] += matching.size();
    if (sharedWire)
	sharedWire->unref();
    if (sharedMesh)
	sharedMesh->unref();

    std::unordered_set<Obol::PartId, std::hash<Obol::PartId>> releasedParts;
    for (const Obol::PartId &oldPart : oldParts) {
	const std::unordered_map<Obol::PartId, size_t,
	    std::hash<Obol::PartId>>::iterator count =
	    this->d->compactIndex->partReferenceCounts.find(oldPart);
	if (count != this->d->compactIndex->partReferenceCounts.end() &&
	    count->second == 0) {
	    releasedParts.insert(oldPart);
	    this->d->compactIndex->partReferenceCounts.erase(count);
	}
    }
    this->d->compactIndex->parts.erase(
	std::remove_if(this->d->compactIndex->parts.begin(),
	    this->d->compactIndex->parts.end(),
	    [&](const BObolCompactPartReference &part) {
		return releasedParts.find(part.part) != releasedParts.end();
	    }), this->d->compactIndex->parts.end());
    for (auto it = this->d->compactIndex->partIdByKey.begin();
	 it != this->d->compactIndex->partIdByKey.end();) {
	if (releasedParts.find(it->second) != releasedParts.end())
	    it = this->d->compactIndex->partIdByKey.erase(it);
	else
	    ++it;
    }
    for (auto it = this->d->compactIndex->partIdByGeometry.begin();
	 it != this->d->compactIndex->partIdByGeometry.end();) {
	if (releasedParts.find(it->second) != releasedParts.end())
	    it = this->d->compactIndex->partIdByGeometry.erase(it);
	else
	    ++it;
    }

    SbBox3f bounds;
    if (compact_index_source_bounds(*this->d->compactIndex, bounds))
	(void)this->setSourceBoundsState(TRUE, bounds.getMin(), bounds.getMax());
    else
	this->clearSourceBounds();
    this->detachFieldSensors();
    this->sourceRevision = revision;
    mark_source_realized_current(this);
    this->attachFieldSensors();
    this->markCompiledAssemblyDirty();
    this->markCadBatchDirty();
    this->markDisplayMeshLodDirty();
    return static_cast<int>(matching.size());
}

static SbMatrix
compact_entry_local_to_world(const BObolCompactInstanceEntry &entry,
			     const SbMatrix &parentToWorld)
{
    SbMatrix matrix = entry.localToSource;
    matrix.multRight(parentToWorld);
    return matrix;
}

template <typename Callback>
static void
compact_for_each_point(const BObolCompactInstanceEntry &entry,
	Callback callback)
{
    if (!entry.geometry || !entry.geometry->points)
	return;
    const Obol::PointRep &points = *entry.geometry->points;
    for (size_t i = 0; i < points.positions.size(); i++) {
	const int primitiveIndex = i < points.pointIds.size() ?
	    static_cast<int>(points.pointIds[i]) : static_cast<int>(i);
	const SbBool colorValid = i < points.colorValid.size() &&
	    points.colorValid[i] && i < points.colors.size();
	const SbBool scaleValid = i < points.scaleValid.size() &&
	    points.scaleValid[i] && i < points.scales.size();
	const SbBool normalValid = i < points.normalValid.size() &&
	    points.normalValid[i] && i < points.normals.size();
	callback(primitiveIndex, points.positions[i], colorValid,
	    colorValid ? points.colors[i] : SbColor(1.0f, 1.0f, 1.0f),
	    scaleValid, scaleValid ? points.scales[i] : 0.0f,
	    normalValid, normalValid ? points.normals[i] :
	    SbVec3f(0.0f, 0.0f, 1.0f));
    }
}

template <typename Callback>
static void
compact_for_each_wire_segment(const BObolCompactInstanceEntry &entry,
	Callback callback)
{
    if (!entry.geometry || !entry.geometry->wire)
	return;
    const Obol::WireRep &wire = *entry.geometry->wire;
    int primitiveIndex = 0;
    for (size_t i = 0; i + 1 < wire.segmentPoints.size(); i += 2)
	callback(primitiveIndex++, wire.segmentPoints[i],
	    wire.segmentPoints[i + 1]);
    for (const Obol::WirePolyline &polyline : wire.polylines) {
	for (size_t i = 1; i < polyline.points.size(); i++)
	    callback(primitiveIndex++, polyline.points[i - 1],
		polyline.points[i]);
    }
}

static float
compact_transform_point_scale(const SbMatrix &matrix, float scale)
{
    SbVec3f x;
    SbVec3f y;
    SbVec3f z;
    matrix.multDirMatrix(SbVec3f(scale, 0.0f, 0.0f), x);
    matrix.multDirMatrix(SbVec3f(0.0f, scale, 0.0f), y);
    matrix.multDirMatrix(SbVec3f(0.0f, 0.0f, scale), z);
    return std::max(x.length(), std::max(y.length(), z.length()));
}

template <typename Callback>
static void
compact_for_each_triangle(const BObolCompactInstanceEntry &entry,
	Callback callback)
{
    if (!entry.geometry || !entry.geometry->shaded)
	return;
    const Obol::TriMesh &mesh = *entry.geometry->shaded;
    int primitiveIndex = 0;
    for (size_t i = 0; i + 2 < mesh.indices.size(); i += 3) {
	const uint32_t ia = mesh.indices[i];
	const uint32_t ib = mesh.indices[i + 1];
	const uint32_t ic = mesh.indices[i + 2];
	if (ia >= mesh.positions.size() || ib >= mesh.positions.size() ||
	    ic >= mesh.positions.size())
	    continue;
	callback(primitiveIndex++, mesh.positions[ia], mesh.positions[ib],
	    mesh.positions[ic], static_cast<int>(ia), static_cast<int>(ib),
	    static_cast<int>(ic));
    }
}

static SbVec3f
compact_closest_point_on_segment(const SbVec3f &query,
				 const SbVec3f &a,
				 const SbVec3f &b)
{
    SbVec3f ab = b - a;
    float len2 = ab.sqrLength();
    if (len2 <= 0.0f)
	return a;

    float t = (query - a).dot(ab) / len2;
    if (t < 0.0f)
	t = 0.0f;
    if (t > 1.0f)
	t = 1.0f;
    return a + ab * t;
}

static SbVec3f
compact_closest_point_on_triangle(const SbVec3f &query,
				  const SbVec3f &a,
				  const SbVec3f &b,
				  const SbVec3f &c)
{
    SbVec3f ab = b - a;
    SbVec3f ac = c - a;
    SbVec3f ap = query - a;
    float d1 = ab.dot(ap);
    float d2 = ac.dot(ap);
    if (d1 <= 0.0f && d2 <= 0.0f)
	return a;

    SbVec3f bp = query - b;
    float d3 = ab.dot(bp);
    float d4 = ac.dot(bp);
    if (d3 >= 0.0f && d4 <= d3)
	return b;

    float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
	float v = d1 / (d1 - d3);
	return a + ab * v;
    }

    SbVec3f cp = query - c;
    float d5 = ab.dot(cp);
    float d6 = ac.dot(cp);
    if (d6 >= 0.0f && d5 <= d6)
	return c;

    float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
	float w = d2 / (d2 - d6);
	return a + ac * w;
    }

    float va = d3 * d6 - d5 * d4;
    if (va <= 0.0f && (d4 - d3) >= 0.0f &&
	(d5 - d6) >= 0.0f) {
	float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
	return b + (c - b) * w;
    }

    float denom = 1.0f / (va + vb + vc);
    float v = vb * denom;
    float w = vc * denom;
    return a + ab * v + ac * w;
}

struct compact_measure_segment_record {
    SbString path;
    SbString editIntentId;
    SbString editIntentRole;
    int primitiveIndex;
    SbVec3f a;
    SbVec3f b;
};

struct compact_measure_angle_record {
    size_t segmentA;
    size_t segmentB;
    SbVec3f shared;
    float degrees;
};

static const float COMPACT_MEASURE_ANGLE_VERTEX_TOLERANCE = 1.0e-5f;

struct compact_measure_endpoint_cell {
    long long x;
    long long y;
    long long z;

    bool operator<(const compact_measure_endpoint_cell &other) const
    {
	if (this->x != other.x)
	    return this->x < other.x;
	if (this->y != other.y)
	    return this->y < other.y;
	return this->z < other.z;
    }
};

static compact_measure_endpoint_cell
compact_measure_make_endpoint_cell(long long x, long long y, long long z)
{
    compact_measure_endpoint_cell cell;
    cell.x = x;
    cell.y = y;
    cell.z = z;
    return cell;
}

static long long
compact_measure_endpoint_cell_coord(float value)
{
    return static_cast<long long>(
	floor(static_cast<double>(value) /
	      static_cast<double>(COMPACT_MEASURE_ANGLE_VERTEX_TOLERANCE)));
}

static compact_measure_endpoint_cell
compact_measure_endpoint_cell_for_point(const SbVec3f &point)
{
    return compact_measure_make_endpoint_cell(
	compact_measure_endpoint_cell_coord(point[0]),
	compact_measure_endpoint_cell_coord(point[1]),
	compact_measure_endpoint_cell_coord(point[2]));
}

static void
compact_measure_collect_angle_endpoint_candidates(
    const std::map<compact_measure_endpoint_cell,
		   std::vector<size_t> > &endpointMap,
    const SbVec3f &point,
    std::vector<size_t> &candidates)
{
    compact_measure_endpoint_cell center =
	compact_measure_endpoint_cell_for_point(point);
    for (int dx = -1; dx <= 1; dx++) {
	for (int dy = -1; dy <= 1; dy++) {
	    for (int dz = -1; dz <= 1; dz++) {
		compact_measure_endpoint_cell key =
		    compact_measure_make_endpoint_cell(center.x + dx,
						       center.y + dy,
						       center.z + dz);
		std::map<compact_measure_endpoint_cell,
		    std::vector<size_t> >::const_iterator it =
			endpointMap.find(key);
		if (it == endpointMap.end())
		    continue;
		candidates.insert(candidates.end(), it->second.begin(),
				  it->second.end());
	    }
	}
    }
}

static void
compact_measure_add_angle_endpoint(
    std::map<compact_measure_endpoint_cell,
	     std::vector<size_t> > &endpointMap,
    const SbVec3f &point,
    size_t segmentIndex)
{
    endpointMap[compact_measure_endpoint_cell_for_point(point)].push_back(
	segmentIndex);
}

static float
compact_measure_clamp_float(float value, float minValue, float maxValue)
{
    if (value < minValue)
	return minValue;
    if (value > maxValue)
	return maxValue;
    return value;
}

static SbBool
compact_measure_same_point(const SbVec3f &a, const SbVec3f &b)
{
    return (a - b).length() <= COMPACT_MEASURE_ANGLE_VERTEX_TOLERANCE;
}

static SbBool
compact_measure_shared_segment_vertex(
    const compact_measure_segment_record &sa,
    const compact_measure_segment_record &sb,
    SbVec3f &shared,
    SbVec3f &otherA,
    SbVec3f &otherB)
{
    if (compact_measure_same_point(sa.a, sb.a)) {
	shared = sa.a;
	otherA = sa.b;
	otherB = sb.b;
	return TRUE;
    }
    if (compact_measure_same_point(sa.a, sb.b)) {
	shared = sa.a;
	otherA = sa.b;
	otherB = sb.a;
	return TRUE;
    }
    if (compact_measure_same_point(sa.b, sb.a)) {
	shared = sa.b;
	otherA = sa.a;
	otherB = sb.b;
	return TRUE;
    }
    if (compact_measure_same_point(sa.b, sb.b)) {
	shared = sa.b;
	otherA = sa.a;
	otherB = sb.a;
	return TRUE;
    }
    return FALSE;
}

static SbBool
compact_measure_segment_angle_degrees(const SbVec3f &shared,
				      const SbVec3f &otherA,
				      const SbVec3f &otherB,
				      float &angleDegrees)
{
    SbVec3f va = otherA - shared;
    SbVec3f vb = otherB - shared;
    float lenA = va.length();
    float lenB = vb.length();
    if (lenA <= 0.0f || lenB <= 0.0f)
	return FALSE;

    float cosAngle = compact_measure_clamp_float(
	va.dot(vb) / (lenA * lenB), -1.0f, 1.0f);
    angleDegrees = acosf(cosAngle) * (180.0f / 3.14159265358979323846f);
    return TRUE;
}

static void
compact_measure_collect_connected_angles(
    const std::vector<compact_measure_segment_record> &measuredSegments,
    std::vector<compact_measure_angle_record> &angleRecords)
{
    if (measuredSegments.size() < 2)
	return;

    std::map<compact_measure_endpoint_cell, std::vector<size_t> > endpointMap;
    std::vector<std::pair<size_t, size_t> > connectedPairs;
    std::vector<size_t> candidates;

    for (size_t i = 0; i < measuredSegments.size(); i++) {
	candidates.clear();
	compact_measure_collect_angle_endpoint_candidates(endpointMap,
	    measuredSegments[i].a, candidates);
	compact_measure_collect_angle_endpoint_candidates(endpointMap,
	    measuredSegments[i].b, candidates);
	std::sort(candidates.begin(), candidates.end());
	candidates.erase(std::unique(candidates.begin(), candidates.end()),
			 candidates.end());
	for (size_t j = 0; j < candidates.size(); j++)
	    connectedPairs.push_back(std::make_pair(candidates[j], i));
	compact_measure_add_angle_endpoint(endpointMap, measuredSegments[i].a,
	    i);
	compact_measure_add_angle_endpoint(endpointMap, measuredSegments[i].b,
	    i);
    }

    std::sort(connectedPairs.begin(), connectedPairs.end());
    connectedPairs.erase(std::unique(connectedPairs.begin(),
				     connectedPairs.end()),
			 connectedPairs.end());
    for (size_t i = 0; i < connectedPairs.size(); i++) {
	size_t segmentA = connectedPairs[i].first;
	size_t segmentB = connectedPairs[i].second;
	SbVec3f shared;
	SbVec3f otherA;
	SbVec3f otherB;
	float degrees = 0.0f;
	if (!compact_measure_shared_segment_vertex(measuredSegments[segmentA],
		measuredSegments[segmentB], shared, otherA, otherB))
	    continue;
	if (!compact_measure_segment_angle_degrees(shared, otherA, otherB,
		degrees))
	    continue;
	compact_measure_angle_record record;
	record.segmentA = segmentA;
	record.segmentB = segmentB;
	record.shared = shared;
	record.degrees = degrees;
	angleRecords.push_back(record);
    }
}

int
SoBRLDatabaseSource::exportCompactInstances(SoBRLExportAction *action,
	const SbMatrix &parentToWorld)
{
    if (!action || !this->hasCompactInstanceIndex())
	return 0;
    if (!this->visible.getValue())
	return 1;

    for (size_t i = 0; i < this->d->compactIndex->entries.size(); i++) {
	const BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[i];
	if (!entry.visible)
	    continue;

	const SbMatrix localToWorld =
	    compact_entry_local_to_world(entry, parentToWorld);
	const BObolRealizedShapeSummary &summary = entry.shapeSummary;
	if (entry.lodBacked &&
	    action->getGeometryPolicy() == SoBRLExportAction::FULL_DETAIL &&
	    entry.sourceMeshRequestValid) {
	    BObolSourceMeshRequest request = entry.sourceMeshRequest;
	    request.localToWorld = localToWorld;
	    action->sourceBackedFullDetailRequests.push_back(request);
	    action->skippedLodDisplayMeshCount++;
	    continue;
	}
	if (entry.pointGeometry) {
	    compact_for_each_point(entry,
		[&](int primitiveIndex, const SbVec3f &localPoint,
		    SbBool pointColorValid, const SbColor &pointColor,
		    SbBool pointScaleValid, float pointScale,
		    SbBool pointNormalValid, const SbVec3f &pointNormal) {
		    SbVec3f worldPoint;
		    SbVec3f worldNormal = pointNormal;
		    localToWorld.multVecMatrix(localPoint, worldPoint);
		    if (pointNormalValid) {
			const SbMatrix normalMatrix =
			    localToWorld.inverse().transpose();
			normalMatrix.multDirMatrix(pointNormal, worldNormal);
			worldNormal.normalize();
		    }
		action->appendPoint(summary.path, summary.sourceName,
			summary.sourceType, summary.sourceId, summary.regionId,
			summary.airCode, summary.materialId, summary.los,
			summary.materialColorValid, summary.materialColor,
			summary.materialShader, primitiveIndex, entry.selected,
			entry.highlighted, summary.ghosted, summary.hiddenLine,
			summary.editEmphasis, summary.editIntentId,
			summary.editIntentRole, summary.lodPolicy,
			summary.colorOverride, summary.color, pointColorValid,
			pointColor, pointScaleValid,
			pointScaleValid ? compact_transform_point_scale(
			    localToWorld, pointScale) : 0.0f,
		    pointNormalValid, worldNormal, worldPoint);
		action->applyLastPointMetadata(summary);
		});
	}
	if (entry.wireGeometry) {
	    compact_for_each_wire_segment(entry,
		[&](int segmentIndex, const SbVec3f &localA,
		    const SbVec3f &localB) {
		SbVec3f worldA;
		SbVec3f worldB;
		localToWorld.multVecMatrix(localA, worldA);
		localToWorld.multVecMatrix(localB, worldB);
		action->appendLine(summary.path, summary.sourceName,
		    summary.sourceType, summary.sourceId, summary.regionId,
		    summary.airCode, summary.materialId, summary.los,
		    summary.materialColorValid, summary.materialColor,
		    summary.materialShader, segmentIndex, entry.selected,
		    entry.highlighted, summary.ghosted, summary.hiddenLine,
		    summary.editEmphasis, summary.lineStyle,
		    summary.lineWidth, summary.editIntentId,
		    summary.editIntentRole, summary.lodPolicy,
		    summary.colorOverride, summary.color,
		    worldA, worldB);
		action->applyLastLineMetadata(summary);
	    });
	    continue;
	}

	if (entry.meshGeometry) {
	    compact_for_each_triangle(entry,
		[&](int triangleIndex, const SbVec3f &a, const SbVec3f &b,
		    const SbVec3f &c, int ia, int ib, int ic) {
		SbVec3f worldA;
		SbVec3f worldB;
		SbVec3f worldC;
		localToWorld.multVecMatrix(a, worldA);
		localToWorld.multVecMatrix(b, worldB);
		localToWorld.multVecMatrix(c, worldC);
		action->appendTriangle(summary.path, summary.sourceName,
		    summary.sourceType, summary.sourceId, summary.regionId,
		    summary.airCode, summary.materialId, summary.los,
		    summary.materialColorValid, summary.materialColor,
		    summary.materialShader, triangleIndex, ia, ib, ic,
		    entry.selected, entry.highlighted, summary.ghosted,
		    summary.hiddenLine, summary.editEmphasis,
		    summary.editIntentId, summary.editIntentRole,
		    summary.lodPolicy, summary.lodAvailable,
		    summary.lodActiveLevel, summary.lodFaceCount,
		    summary.lodPointCount, summary.lodOriginalPointCount,
		    summary.lodNormalCount, summary.lodHasSnappedPoints,
		    summary.lodHasNormals, summary.lodBoundsMin,
		    summary.lodBoundsMax, summary.colorOverride, summary.color,
		    worldA, worldB, worldC);
		action->applyLastTriangleMetadata(summary);
	    });
	}
    }
    return 1;
}

int
SoBRLDatabaseSource::measureCompactInstances(SoBRLMeasureAction *action,
	const SbMatrix &parentToWorld)
{
    if (!action || !this->hasCompactInstanceIndex())
	return 0;
    if (!this->visible.getValue())
	return 1;

    for (size_t i = 0; i < this->d->compactIndex->entries.size(); i++) {
	const BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[i];
	if (!entry.visible ||
	    !action->selectionAllows(entry.selected) ||
	    !action->highlightAllows(entry.highlighted))
	    continue;

	const SbMatrix localToWorld =
	    compact_entry_local_to_world(entry, parentToWorld);
	SbBool measuredShape = FALSE;
	if (entry.wireGeometry) {
	    std::vector<compact_measure_segment_record> measuredSegments;
	    compact_for_each_wire_segment(entry,
		[&](int segmentIndex, const SbVec3f &localA,
		    const SbVec3f &localB) {
		SbVec3f pointA =
		    action->pointForCoordinateSpace(localToWorld, localA);
		SbVec3f pointB =
		    action->pointForCoordinateSpace(localToWorld, localB);
		action->measureSegment(entry.semantic.path,
		    entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, segmentIndex,
		    pointA, pointB);
		if (action->angleComputationEnabled) {
		    compact_measure_segment_record record;
		    record.path = entry.semantic.path;
		    record.editIntentId = entry.semantic.editIntentId;
		    record.editIntentRole = entry.semantic.editIntentRole;
		    record.primitiveIndex = segmentIndex;
		    record.a = pointA;
		    record.b = pointB;
		    measuredSegments.push_back(record);
		}
		measuredShape = TRUE;
	    });
	    if (action->angleComputationEnabled) {
		std::vector<compact_measure_angle_record> angleRecords;
		compact_measure_collect_connected_angles(measuredSegments,
		    angleRecords);
		for (size_t j = 0; j < angleRecords.size(); j++) {
		    const compact_measure_angle_record &angle =
			angleRecords[j];
		    const compact_measure_segment_record &segmentA =
			measuredSegments[angle.segmentA];
		    const compact_measure_segment_record &segmentB =
			measuredSegments[angle.segmentB];
		    action->considerAngle(segmentA.path, segmentA.editIntentId,
			segmentA.editIntentRole, segmentA.primitiveIndex,
			segmentB.primitiveIndex, angle.shared, angle.degrees);
		}
	    }
	} else if (entry.meshGeometry) {
	    compact_for_each_triangle(entry,
		[&](int triangleIndex, const SbVec3f &a, const SbVec3f &b,
		    const SbVec3f &c, int ia, int ib, int ic) {
		SbVec3f pointA =
		    action->pointForCoordinateSpace(localToWorld, a);
		SbVec3f pointB =
		    action->pointForCoordinateSpace(localToWorld, b);
		SbVec3f pointC =
		    action->pointForCoordinateSpace(localToWorld, c);
		action->measureTriangle(entry.semantic.path,
		    entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, pointA, pointB,
		    pointC, ia, ib, ic);
		measuredShape = TRUE;
	    });
	}
	if (measuredShape)
	    action->shapeCount++;
    }
    return 1;
}

int
SoBRLDatabaseSource::snapCompactInstances(SoBRLSnapAction *action,
	const SbMatrix &parentToWorld)
{
    if (!action || !this->hasCompactInstanceIndex())
	return 0;
    if (!this->visible.getValue())
	return 1;

    const SbVec3f query = action->queryPoint;
    for (size_t i = 0; i < this->d->compactIndex->entries.size(); i++) {
	const BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[i];
	if (!entry.visible || !entry.selectable ||
	    !action->selectionAllows(entry.selected))
	    continue;

	const SbMatrix localToWorld =
	    compact_entry_local_to_world(entry, parentToWorld);
	SbBox3f centerBox;
	centerBox.makeEmpty();
	if (entry.pointGeometry) {
	    compact_for_each_point(entry,
		[&](int primitiveIndex, const SbVec3f &localPoint,
		    SbBool, const SbColor &, SbBool, float,
		    SbBool, const SbVec3f &) {
		    const SbVec3f point = action->pointForCoordinateSpace(
			localToWorld, localPoint);
		    centerBox.extendBy(point);
		    action->consider(SoBRLSnapAction::ENDPOINT,
			entry.semantic.path, entry.semantic.editIntentId,
			entry.semantic.editIntentRole, primitiveIndex, query,
			point);
		});
	}
	if (entry.wireGeometry) {
	    compact_for_each_wire_segment(entry,
		[&](int segmentIndex, const SbVec3f &localA,
		    const SbVec3f &localB) {
		SbVec3f pointA =
		    action->pointForCoordinateSpace(localToWorld, localA);
		SbVec3f pointB =
		    action->pointForCoordinateSpace(localToWorld, localB);
		centerBox.extendBy(pointA);
		centerBox.extendBy(pointB);
		action->consider(SoBRLSnapAction::ENDPOINT,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, segmentIndex, query,
		    pointA);
		action->consider(SoBRLSnapAction::ENDPOINT,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, segmentIndex, query,
		    pointB);
		action->consider(SoBRLSnapAction::LINE_NEAREST,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, segmentIndex, query,
		    compact_closest_point_on_segment(query, pointA, pointB));
		action->consider(SoBRLSnapAction::MIDPOINT,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, segmentIndex, query,
		    (pointA + pointB) * 0.5f);
	    });
	} else if (entry.meshGeometry) {
	    compact_for_each_triangle(entry,
		[&](int triangleIndex, const SbVec3f &a, const SbVec3f &b,
		    const SbVec3f &c, int ia, int ib, int ic) {
		SbVec3f pointA =
		    action->pointForCoordinateSpace(localToWorld, a);
		SbVec3f pointB =
		    action->pointForCoordinateSpace(localToWorld, b);
		SbVec3f pointC =
		    action->pointForCoordinateSpace(localToWorld, c);
		centerBox.extendBy(pointA);
		centerBox.extendBy(pointB);
		centerBox.extendBy(pointC);
		action->consider(SoBRLSnapAction::VERTEX,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, query, pointA, ia);
		action->consider(SoBRLSnapAction::VERTEX,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, query, pointB, ib);
		action->consider(SoBRLSnapAction::VERTEX,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, query, pointC, ic);
		action->consider(SoBRLSnapAction::FACE_NEAREST,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, query,
		    compact_closest_point_on_triangle(query, pointA, pointB,
					      pointC));
		action->consider(SoBRLSnapAction::EDGE_NEAREST,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, query,
		    compact_closest_point_on_segment(query, pointA, pointB),
		    -1, 0, ia, ib);
		action->consider(SoBRLSnapAction::EDGE_NEAREST,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, query,
		    compact_closest_point_on_segment(query, pointB, pointC),
		    -1, 1, ib, ic);
		action->consider(SoBRLSnapAction::EDGE_NEAREST,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, query,
		    compact_closest_point_on_segment(query, pointC, pointA),
		    -1, 2, ic, ia);
	    });
	}
	if (!centerBox.isEmpty()) {
	    const SoBRLCadAssembly::InstanceSemantic &semantic =
		entry.semantic;
	    action->consider(SoBRLSnapAction::CENTER, semantic.path,
		semantic.editIntentId, semantic.editIntentRole, -1, query,
		centerBox.getCenter());
	}
    }
    return 1;
}

int
SoBRLDatabaseSource::prepareCompiledAssembly(void)
{
    return this->syncCompiledAssembly();
}

SbBool
SoBRLDatabaseSource::hasCompiledAssembly(void) const
{
    return (this->d->compiledAssembly && this->d->compiledAssemblyActive) ?
	   TRUE : FALSE;
}

int
SoBRLDatabaseSource::getCompiledAssemblyPartCount(void) const
{
    if (!this->d->compiledAssembly || !this->d->compiledAssemblyActive)
	return 0;
    return static_cast<int>(this->d->compiledAssembly->partCount());
}

int
SoBRLDatabaseSource::getCompiledAssemblyInstanceCount(void) const
{
    if (!this->d->compiledAssembly || !this->d->compiledAssemblyActive)
	return 0;
    return static_cast<int>(this->d->compiledAssembly->instanceCount());
}

static int
count_material_objects_in_node(SoNode *node)
{
    if (!node)
	return 0;

    if (node->isOfType(SoBRLMaterialObject::getClassTypeId()))
	return 1;

    int ret = 0;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    ret += count_material_objects_in_node(group->getChild(i));
    }

    return ret;
}

int
SoBRLDatabaseSource::getRealizedMaterialObjectCount(void) const
{
    int ret = 0;
    for (int i = 0; i < this->getNumChildren(); i++)
	ret += count_material_objects_in_node(this->getChild(i));
    return ret;
}

static void
realized_material_summary_owner(const SoBRLDatabaseSource *source,
				int sourceIndex, BObolRealizedMaterialSummary &summary);

static void
realized_material_summary(const SoBRLMaterialObject *object,
			  BObolRealizedMaterialSummary &summary)
{
    summary = BObolRealizedMaterialSummary();
    if (!object)
	return;

    summary.valid = TRUE;
    summary.sourcePath = object->sourcePath.getValue();
    summary.sourceName = object->sourceName.getValue();
    summary.sourceType = object->sourceType.getValue();
    summary.sourceId = object->sourceId.getValue();
    summary.materialName = object->materialName.getValue();
    summary.parentName = object->parentName.getValue();
    summary.materialSource = object->materialSource.getValue();
    summary.propertyCount = object->getPropertyCount();
}

int
SoBRLDatabaseSource::getRealizedMaterialSummaryCount(void) const
{
    return this->getRealizedMaterialObjectCount();
}

SbBool
SoBRLDatabaseSource::getRealizedMaterialSummary(int index,
	BObolRealizedMaterialSummary &summary) const
{
    summary = BObolRealizedMaterialSummary();
    if (index < 0)
	return FALSE;

    SoBRLMaterialObject *object = this->getRealizedMaterialObject(index);
    if (!object)
	return FALSE;

    realized_material_summary(object, summary);
    realized_material_summary_owner(this, -1, summary);
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getRealizedMaterialProperty(int materialIndex,
	int propertyIndex, SbString &groupOut, SbString &nameOut,
	SbString &valueOut) const
{
    if (materialIndex < 0 || propertyIndex < 0)
	return FALSE;

    SoBRLMaterialObject *object =
	this->getRealizedMaterialObject(materialIndex);
    if (!object)
	return FALSE;

    return object->getProperty(propertyIndex, groupOut, nameOut, valueOut);
}

static void
realized_summary_owner_instance_from_node(const SoNode *node,
	SbString &ownerSourceInstanceKey)
{
    if (!node || ownerSourceInstanceKey.getLength() > 0)
	return;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	ownerSourceInstanceKey = source_effective_instance_key(source);
	return;
    }
    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	ownerSourceInstanceKey =
	    static_cast<const SoBRLVListShape *>(node)->
	    ownerSourceInstanceKey.getValue();
	return;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	ownerSourceInstanceKey =
	    static_cast<const SoBRLMeshShape *>(node)->
	    ownerSourceInstanceKey.getValue();
	return;
    }
}

static void
realized_tree_summary_fill(const SoNode *node, int depth, SbBool hasParent,
			   int ownerSourceIndex, const SbString &ownerSourcePath,
			   BObolSceneTreeSummary &summary)
{
    summary = BObolSceneTreeSummary();
    if (!node)
	return;

    summary.valid = TRUE;
    summary.hasParent = hasParent;
    summary.drawTreeDepth = depth;
    summary.ownerSourceIndex = ownerSourceIndex;
    summary.ownerSourcePath = ownerSourcePath;
    realized_summary_owner_instance_from_node(node,
	    summary.ownerSourceInstanceKey);
    summary.isGroup =
	node->isOfType(SoGroup::getClassTypeId()) ? TRUE : FALSE;
    summary.isDatabaseSource =
	node->isOfType(SoBRLDatabaseSource::getClassTypeId()) ? TRUE : FALSE;
    summary.isMaterialObject =
	node->isOfType(SoBRLMaterialObject::getClassTypeId()) ? TRUE : FALSE;
    summary.isShape =
	(node->isOfType(SoBRLVListShape::getClassTypeId()) ||
	 node->isOfType(SoBRLMeshShape::getClassTypeId())) ? TRUE : FALSE;

    if (summary.isGroup)
	summary.childCount = static_cast<const SoGroup *>(node)->getNumChildren();

    if (summary.isDatabaseSource) {
	const SoBRLDatabaseSource *sourceNode =
	    static_cast<const SoBRLDatabaseSource *>(node);
	summary.nodeKind = BObolSceneTreeSummary::NODE_DATABASE_SOURCE;
	summary.path = sourceNode->path.getValue();
	summary.displayName = sourceNode->displayName.getValue();
	if (summary.ownerSourcePath.getLength() == 0)
	    summary.ownerSourcePath = sourceNode->path.getValue();
	if (summary.ownerSourceInstanceKey.getLength() == 0)
	    summary.ownerSourceInstanceKey =
		source_effective_instance_key(sourceNode);
	return;
    }

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	const SoBRLVListShape *shape =
	    static_cast<const SoBRLVListShape *>(node);
	summary.nodeKind = BObolSceneTreeSummary::NODE_VLIST_SHAPE;
	summary.path = shape->sourcePath.getValue();
	summary.sourceName = shape->sourceName.getValue();
	summary.sourceType = shape->sourceType.getValue();
	summary.sourceId = shape->sourceId.getValue();
	summary.displayName = shape->displayName.getValue();
	summary.geometryName = shape->geometryName.getValue();
	return;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	const SoBRLMeshShape *shape =
	    static_cast<const SoBRLMeshShape *>(node);
	summary.nodeKind = BObolSceneTreeSummary::NODE_MESH_SHAPE;
	summary.path = shape->sourcePath.getValue();
	summary.sourceName = shape->sourceName.getValue();
	summary.sourceType = shape->sourceType.getValue();
	summary.sourceId = shape->sourceId.getValue();
	summary.displayName = shape->displayName.getValue();
	summary.geometryName = shape->geometryName.getValue();
	return;
    }

    if (summary.isMaterialObject) {
	const SoBRLMaterialObject *object =
	    static_cast<const SoBRLMaterialObject *>(node);
	summary.nodeKind = BObolSceneTreeSummary::NODE_MATERIAL_OBJECT;
	summary.path = object->sourcePath.getValue();
	summary.sourceName = object->sourceName.getValue();
	summary.sourceType = object->sourceType.getValue();
	summary.sourceId = object->sourceId.getValue();
	summary.displayName = object->materialName.getValue();
	return;
    }

    summary.nodeKind = summary.isGroup ?
		       BObolSceneTreeSummary::NODE_GROUP :
		       BObolSceneTreeSummary::NODE_OTHER;
}

static int
realized_tree_summary_node_count(const SoNode *node)
{
    if (!node)
	return 0;
    if (node_is_source_placement_transform(node))
	return 0;

    int count = 1;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    count += realized_tree_summary_node_count(group->getChild(i));
    }
    return count;
}

static SbBool
find_realized_tree_summary_in_node(const SoNode *node, int &index, int depth,
				   SbBool hasParent, int ownerSourceIndex,
				   const SbString &ownerSourcePath, BObolSceneTreeSummary &summary)
{
    if (!node)
	return FALSE;
    if (node_is_source_placement_transform(node))
	return FALSE;

    if (index == 0) {
	realized_tree_summary_fill(node, depth, hasParent, ownerSourceIndex,
				   ownerSourcePath, summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    if (find_realized_tree_summary_in_node(group->getChild(i),
						   index, depth + 1, TRUE, ownerSourceIndex,
						   ownerSourcePath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLDatabaseSource::getRealizedTreeSummaryCount(void) const
{
    return realized_tree_summary_node_count(this);
}

SbBool
SoBRLDatabaseSource::getRealizedTreeSummary(int index,
	BObolSceneTreeSummary &summary) const
{
    summary = BObolSceneTreeSummary();
    if (index < 0)
	return FALSE;

    const SbString ownerPath = this->path.getValue();
    const SbBool ret = find_realized_tree_summary_in_node(this, index, 0,
		       FALSE, -1, ownerPath, summary);
    if (ret && summary.ownerSourceInstanceKey.getLength() == 0)
	summary.ownerSourceInstanceKey = source_effective_instance_key(this);
    return ret;
}

static void
realized_display_summary_fill_common(BObolSceneDisplaySummary &summary,
				     int nodeKind, int ownerSourceIndex, const SbString &ownerSourcePath,
				     const SbString &nodePath)
{
    summary = BObolSceneDisplaySummary();
    summary.valid = TRUE;
    summary.nodeKind = nodeKind;
    summary.ownerSourceIndex = ownerSourceIndex;
    summary.ownerSourcePath = ownerSourcePath;
    summary.path = nodePath;
}

template <typename ShapeT>
static void
realized_display_summary_fill_shape(const ShapeT *shape,
				    int nodeKind, int ownerSourceIndex, const SbString &ownerSourcePath,
				    BObolSceneDisplaySummary &summary)
{
    realized_display_summary_fill_common(summary, nodeKind, ownerSourceIndex,
					 ownerSourcePath, shape ? shape->sourcePath.getValue() : SbString(""));
    if (!shape)
	return;

    summary.ownerSourceInstanceKey = shape->ownerSourceInstanceKey.getValue();
    summary.hasDrawIntent = shape->sourcePath.getValue().getLength() > 0;
    summary.intentPath = shape->sourcePath.getValue();
    summary.intentDrawMode = shape->drawMode.getValue();
    summary.visible = shape->visible.getValue();
    summary.highlighted = shape->highlighted.getValue();
    summary.lineStyle = shape->lineStyle.getValue();
    summary.lineWidth = shape->lineWidth.getValue();
    summary.transparency = shape->transparency.getValue();
    summary.drawMode = shape->drawMode.getValue();
    summary.materialValid = TRUE;
    summary.materialRevision = shape->materialRevision.getValue();
    if (shape->materialColorValid.getValue())
	summary.materialColor = shape->materialColor.getValue();
    else if (shape->colorOverride.getValue())
	summary.materialColor = shape->color.getValue();
    summary.drawMatrixValid = shape->drawMatrixValid.getValue();
    summary.drawMatrix = shape->drawMatrix.getValue();
    summary.drawCenterValid = shape->drawCenterValid.getValue();
    summary.drawCenter = shape->drawCenter.getValue();
    summary.drawSizeValid = shape->drawSizeValid.getValue();
    summary.drawSize = shape->drawSize.getValue();
}

static void
realized_display_summary_fill(const SoNode *node, int ownerSourceIndex,
			      const SbString &ownerSourcePath, BObolSceneDisplaySummary &summary)
{
    summary = BObolSceneDisplaySummary();
    if (!node)
	return;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	realized_display_summary_fill_common(summary,
					     BObolSceneTreeSummary::NODE_DATABASE_SOURCE,
					     ownerSourceIndex, ownerSourcePath, source->path.getValue());
	summary.ownerSourceInstanceKey = source_effective_instance_key(source);
	summary.isDatabaseSource = TRUE;
	summary.hasDrawIntent = source->path.getValue().getLength() > 0;
	summary.intentPath = source->path.getValue();
	summary.intentDrawMode = source->drawMode.getValue();
	summary.visible = source->visible.getValue();
	summary.selected = source->selected.getValue();
	summary.highlighted = source->highlighted.getValue();
	summary.lineStyle = source->lineStyle.getValue();
	summary.lineWidth = source->lineWidth.getValue();
	summary.transparency = source->transparency.getValue();
	summary.drawMode = source->drawMode.getValue();
	summary.materialValid = TRUE;
	summary.materialRevision = source->materialRevision.getValue();
	if (source->materialColorValid.getValue())
	    summary.materialColor = source->materialColor.getValue();
	else if (source->colorOverride.getValue())
	    summary.materialColor = source->color.getValue();
	return;
    }

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	realized_display_summary_fill_shape(
	    static_cast<const SoBRLVListShape *>(node),
	    BObolSceneTreeSummary::NODE_VLIST_SHAPE,
	    ownerSourceIndex, ownerSourcePath, summary);
	return;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	realized_display_summary_fill_shape(
	    static_cast<const SoBRLMeshShape *>(node),
	    BObolSceneTreeSummary::NODE_MESH_SHAPE,
	    ownerSourceIndex, ownerSourcePath, summary);
	return;
    }

    if (node->isOfType(SoBRLMaterialObject::getClassTypeId())) {
	const SoBRLMaterialObject *object =
	    static_cast<const SoBRLMaterialObject *>(node);
	realized_display_summary_fill_common(summary,
					     BObolSceneTreeSummary::NODE_MATERIAL_OBJECT,
					     ownerSourceIndex, ownerSourcePath,
					     object->sourcePath.getValue());
	return;
    }

    const int nodeKind = node->isOfType(SoGroup::getClassTypeId()) ?
			 BObolSceneTreeSummary::NODE_GROUP :
			 BObolSceneTreeSummary::NODE_OTHER;
    realized_display_summary_fill_common(summary, nodeKind, ownerSourceIndex,
					 ownerSourcePath, "");
}

static SbBool
find_realized_display_summary_in_node(const SoNode *node, int &index,
				      int ownerSourceIndex, const SbString &ownerSourcePath,
				      BObolSceneDisplaySummary &summary)
{
    if (!node)
	return FALSE;
    if (node_is_source_placement_transform(node))
	return FALSE;

    if (index == 0) {
	realized_display_summary_fill(node, ownerSourceIndex,
				      ownerSourcePath, summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    if (find_realized_display_summary_in_node(group->getChild(i),
		index, ownerSourceIndex, ownerSourcePath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLDatabaseSource::getRealizedDisplaySummaryCount(void) const
{
    return this->getRealizedTreeSummaryCount();
}

SbBool
SoBRLDatabaseSource::getRealizedDisplaySummary(int index,
	BObolSceneDisplaySummary &summary) const
{
    summary = BObolSceneDisplaySummary();
    if (index < 0)
	return FALSE;

    const SbString ownerPath = this->path.getValue();
    const SbBool ret = find_realized_display_summary_in_node(this, index,
		       -1, ownerPath, summary);
    if (ret && summary.ownerSourceInstanceKey.getLength() == 0)
	summary.ownerSourceInstanceKey = source_effective_instance_key(this);
    return ret;
}

static void
scene_material_summary_from_display(const BObolSceneDisplaySummary &display,
				    BObolSceneMaterialSummary &summary)
{
    summary = BObolSceneMaterialSummary();
    if (!display.valid)
	return;

    summary.valid = TRUE;
    summary.nodeKind = display.nodeKind;
    summary.materialValid =
	(display.nodeKind == BObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	 display.nodeKind == BObolSceneTreeSummary::NODE_MESH_SHAPE) ?
	display.materialValid : FALSE;
    summary.materialRevision = display.materialRevision;
    summary.materialColor = display.materialColor;
    summary.ownerSourceIndex = display.ownerSourceIndex;
    summary.ownerSourcePath = display.ownerSourcePath;
    summary.ownerSourceInstanceKey = display.ownerSourceInstanceKey;
    summary.path = display.path;
}

int
SoBRLDatabaseSource::getRealizedSceneMaterialSummaryCount(void) const
{
    return this->getRealizedDisplaySummaryCount();
}

SbBool
SoBRLDatabaseSource::getRealizedSceneMaterialSummary(int index,
	BObolSceneMaterialSummary &summary) const
{
    summary = BObolSceneMaterialSummary();
    BObolSceneDisplaySummary display;
    if (!this->getRealizedDisplaySummary(index, display))
	return FALSE;

    scene_material_summary_from_display(display, summary);
    return TRUE;
}

static int
realized_bounds_node_kind(const SoNode *node)
{
    if (!node)
	return BObolSceneTreeSummary::NODE_UNKNOWN;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	return BObolSceneTreeSummary::NODE_DATABASE_SOURCE;
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return BObolSceneTreeSummary::NODE_VLIST_SHAPE;
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return BObolSceneTreeSummary::NODE_MESH_SHAPE;
    if (node->isOfType(SoBRLMaterialObject::getClassTypeId()))
	return BObolSceneTreeSummary::NODE_MATERIAL_OBJECT;
    if (node->isOfType(SoGroup::getClassTypeId()))
	return BObolSceneTreeSummary::NODE_GROUP;
    return BObolSceneTreeSummary::NODE_OTHER;
}

static SbString
realized_bounds_node_path(const SoNode *node)
{
    if (!node)
	return "";
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	return static_cast<const SoBRLDatabaseSource *>(node)->path.getValue();
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<const SoBRLVListShape *>(node)->sourcePath.getValue();
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<const SoBRLMeshShape *>(node)->sourcePath.getValue();
    if (node->isOfType(SoBRLMaterialObject::getClassTypeId()))
	return static_cast<const SoBRLMaterialObject *>(node)->sourcePath.getValue();
    return "";
}

static SbBool
realized_bounds_for_vlist_shape(const SoBRLVListShape *shape, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!shape)
	return FALSE;
    const SoBRLVListShape *geom = shape->getGeometrySource();
    for (int i = 0; i < geom->point.getNum(); i++)
	bounds.extendBy(geom->point[i]);
    return geom->point.getNum() > 0;
}

static SbBool
realized_bounds_for_mesh_shape(const SoBRLMeshShape *shape, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!shape)
	return FALSE;
    const SoBRLMeshShape *geom = shape->getGeometrySource();
    for (int i = 0; i < geom->point.getNum(); i++)
	bounds.extendBy(geom->point[i]);
    return geom->point.getNum() > 0;
}

static SbBool
realized_bounds_for_node(const SoNode *node, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!node)
	return FALSE;

    SoGetBoundingBoxAction bboxAction(SbViewportRegion(1, 1));
    bboxAction.apply(const_cast<SoNode *>(node));
    bounds = bboxAction.getBoundingBox();
    if (!bounds.isEmpty())
	return TRUE;

    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return realized_bounds_for_vlist_shape(
		   static_cast<const SoBRLVListShape *>(node), bounds);

    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return realized_bounds_for_mesh_shape(
		   static_cast<const SoBRLMeshShape *>(node), bounds);

    SbBool valid = FALSE;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    SbBox3f childBounds;
	    if (realized_bounds_for_node(group->getChild(i),
					 childBounds)) {
		bounds.extendBy(childBounds);
		valid = TRUE;
	    }
	}
    }

    if (valid)
	return TRUE;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()) &&
	static_cast<const SoBRLDatabaseSource *>(node)->
	getEffectiveSourceBounds(bounds))
	return TRUE;

    return valid;
}

static void
realized_bounds_summary_fill(const SoNode *node, int ownerSourceIndex,
			     const SbString &ownerSourcePath, BObolSceneBoundsSummary &summary)
{
    summary = BObolSceneBoundsSummary();
    if (!node)
	return;

    summary.valid = TRUE;
    summary.nodeKind = realized_bounds_node_kind(node);
    summary.ownerSourceIndex = ownerSourceIndex;
    summary.ownerSourcePath = ownerSourcePath;
    realized_summary_owner_instance_from_node(node,
	    summary.ownerSourceInstanceKey);
    summary.path = realized_bounds_node_path(node);
    summary.boundsValid = realized_bounds_for_node(node, summary.bounds);
}

static SbBool
find_realized_bounds_summary_in_node(const SoNode *node, int &index,
				     int ownerSourceIndex, const SbString &ownerSourcePath,
				     BObolSceneBoundsSummary &summary)
{
    if (!node)
	return FALSE;
    if (node_is_source_placement_transform(node))
	return FALSE;

    if (index == 0) {
	realized_bounds_summary_fill(node, ownerSourceIndex, ownerSourcePath,
				     summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    if (find_realized_bounds_summary_in_node(group->getChild(i),
		index, ownerSourceIndex, ownerSourcePath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLDatabaseSource::getRealizedBoundsSummaryCount(void) const
{
    return this->getRealizedTreeSummaryCount();
}

SbBool
SoBRLDatabaseSource::getRealizedBoundsSummary(int index,
	BObolSceneBoundsSummary &summary) const
{
    summary = BObolSceneBoundsSummary();
    if (index < 0)
	return FALSE;

    const SbString ownerPath = this->path.getValue();
    const SbBool ret = find_realized_bounds_summary_in_node(this, index,
		       -1, ownerPath, summary);
    if (ret && summary.ownerSourceInstanceKey.getLength() == 0)
	summary.ownerSourceInstanceKey = source_effective_instance_key(this);
    return ret;
}

template <typename ShapeT>
static void
realized_shape_summary_common(const ShapeT *shape,
			      BObolRealizedShapeSummary &summary)
{
    summary.path = shape->sourcePath.getValue();
    summary.sourceName = shape->sourceName.getValue();
    summary.sourceType = shape->sourceType.getValue();
    summary.sourceId = shape->sourceId.getValue();
    summary.displayName = shape->displayName.getValue();
    summary.geometryName = shape->geometryName.getValue();
    summary.cacheIdentity = shape->cacheIdentity.getValue();
    summary.sourceIdentity = shape->sourceIdentity.getValue();
    summary.databaseIntent = shape->databaseIntent.getValue();
    summary.overlayIntent = shape->overlayIntent.getValue();
    summary.hudIntent = shape->hudIntent.getValue();
    summary.localSource = shape->localSource.getValue();
    summary.sharedSource = shape->sharedSource.getValue();
    summary.nonDatabaseSource = shape->nonDatabaseSource.getValue();
    summary.drawMode = shape->drawMode.getValue();
    summary.recordRole = shape->recordRole.getValue();
    summary.geometryKind = shape->geometryKind.getValue();
    summary.regionId = shape->regionId.getValue();
    summary.airCode = shape->airCode.getValue();
    summary.materialId = shape->materialId.getValue();
    summary.los = shape->los.getValue();
    summary.materialColorValid = shape->materialColorValid.getValue();
    summary.materialColor = shape->materialColor.getValue();
    summary.materialShader = shape->materialShader.getValue();
    summary.materialRevision = shape->materialRevision.getValue();
    summary.visible = shape->visible.getValue();
    summary.selectable = shape->selectable.getValue();
    summary.selected = shape->selected.getValue();
    summary.highlighted = shape->highlighted.getValue();
    summary.ghosted = shape->ghosted.getValue();
    summary.hiddenLine = shape->hiddenLine.getValue();
    summary.editEmphasis = shape->editEmphasis.getValue();
    summary.lineStyle = shape->lineStyle.getValue();
    summary.lineWidth = shape->lineWidth.getValue();
    summary.transparency = shape->transparency.getValue();
    summary.editIntentId = shape->editIntentId.getValue();
    summary.editIntentRole = shape->editIntentRole.getValue();
    summary.lodPolicy = shape->lodPolicy.getValue();
    summary.colorOverride = shape->colorOverride.getValue();
    summary.color = shape->color.getValue();
}

static void
realized_shape_summary_owner(const SoBRLDatabaseSource *source,
			     int sourceIndex, BObolRealizedShapeSummary &summary)
{
    if (!source)
	return;

    summary.ownerSourceIndex = sourceIndex;
    summary.ownerSourcePath = source->path.getValue();
    summary.ownerSourceInstanceKey = source_effective_instance_key(source);
    summary.ownerDrawMode = source->drawMode.getValue();
    summary.ownerSourceRevision = source->sourceRevision.getValue();
    summary.ownerInputsRevision = source->inputsRevision.getValue();
    summary.ownerViewRevision = source->viewRevision.getValue();
    summary.ownerRealizedRevision = source->realizedRevision.getValue();
    summary.ownerRealizedSourceRevision =
	source->realizedSourceRevision.getValue();
    summary.ownerRealizedInputsRevision =
	source->realizedInputsRevision.getValue();
    summary.ownerRealizedViewRevision =
	source->realizedViewRevision.getValue();
    summary.ownerRealizationStatus = source->realizationStatus.getValue();
    summary.ownerRealizationDiagnostic =
	source->realizationDiagnostic.getValue();
    summary.ownerRealizationIdentity =
	source->realizationIdentity.getValue();
    summary.ownerSourceStale = source->stale.getValue();
    summary.ownerStaleReason = source->staleReason.getValue();
}

static void
realized_material_summary_owner(const SoBRLDatabaseSource *source,
				int sourceIndex, BObolRealizedMaterialSummary &summary)
{
    if (!source)
	return;

    summary.ownerSourceIndex = sourceIndex;
    summary.ownerSourcePath = source->path.getValue();
    summary.ownerSourceInstanceKey = source_effective_instance_key(source);
    summary.ownerDrawMode = source->drawMode.getValue();
    summary.ownerSourceRevision = source->sourceRevision.getValue();
    summary.ownerInputsRevision = source->inputsRevision.getValue();
    summary.ownerViewRevision = source->viewRevision.getValue();
    summary.ownerRealizedRevision = source->realizedRevision.getValue();
    summary.ownerRealizedSourceRevision =
	source->realizedSourceRevision.getValue();
    summary.ownerRealizedInputsRevision =
	source->realizedInputsRevision.getValue();
    summary.ownerRealizedViewRevision =
	source->realizedViewRevision.getValue();
    summary.ownerRealizationStatus = source->realizationStatus.getValue();
    summary.ownerRealizationDiagnostic =
	source->realizationDiagnostic.getValue();
    summary.ownerRealizationIdentity =
	source->realizationIdentity.getValue();
    summary.ownerSourceStale = source->stale.getValue();
    summary.ownerStaleReason = source->staleReason.getValue();
}

static void
realized_shape_summary_bounds(const SoMFVec3f &points,
			      BObolRealizedShapeSummary &summary)
{
    summary.bounds.makeEmpty();
    summary.boundsValid = FALSE;
    for (int i = 0; i < points.getNum(); i++) {
	summary.bounds.extendBy(points[i]);
	summary.boundsValid = TRUE;
    }
}

static void
realized_vlist_shape_summary(const SoBRLVListShape *shape,
			     BObolRealizedShapeSummary &summary)
{
    summary = BObolRealizedShapeSummary();
    if (!shape)
	return;

    summary.valid = TRUE;
    summary.shapeKind = BObolRealizedShapeSummary::SHAPE_VLIST;
    realized_shape_summary_common(shape, summary);
    const SoBRLVListShape *geom = shape->getGeometrySource();
    summary.pointCount = geom->point.getNum();
    summary.commandCount = geom->command.getNum();
    summary.segmentCount = shape->getSegmentCount();
    summary.pointPrimitiveCount = shape->getPointPrimitiveCount();
    realized_shape_summary_bounds(geom->point, summary);
}

static void
realized_mesh_shape_summary(const SoBRLMeshShape *shape,
			    BObolRealizedShapeSummary &summary)
{
    summary = BObolRealizedShapeSummary();
    if (!shape)
	return;

    summary.valid = TRUE;
    summary.shapeKind = BObolRealizedShapeSummary::SHAPE_MESH;
    realized_shape_summary_common(shape, summary);
    const SoBRLMeshShape *geom = shape->getGeometrySource();
    summary.pointCount = geom->point.getNum();
    summary.indexCount = geom->coordIndex.getNum();
    summary.triangleCount = shape->getTriangleCount();
    realized_shape_summary_bounds(geom->point, summary);
    summary.lodAvailable = shape->lodAvailable.getValue();
    summary.lodActiveLevel = shape->lodActiveLevel.getValue();
    summary.lodFaceCount = shape->lodFaceCount.getValue();
    summary.lodPointCount = shape->lodPointCount.getValue();
    summary.lodOriginalPointCount = shape->lodOriginalPointCount.getValue();
    summary.lodNormalCount = shape->lodNormalCount.getValue();
    summary.lodHasSnappedPoints = shape->lodHasSnappedPoints.getValue();
    summary.lodHasNormals = shape->lodHasNormals.getValue();
    summary.lodBoundsMin = shape->lodBoundsMin.getValue();
    summary.lodBoundsMax = shape->lodBoundsMax.getValue();
}

static SbBool
find_realized_shape_summary_in_node(SoNode *node, int &index,
				    BObolRealizedShapeSummary &summary)
{
    if (!node)
	return FALSE;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	if (index == 0) {
	    realized_vlist_shape_summary(
		static_cast<SoBRLVListShape *>(node), summary);
	    return TRUE;
	}
	index--;
	return FALSE;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	if (index == 0) {
	    realized_mesh_shape_summary(
		static_cast<SoBRLMeshShape *>(node), summary);
	    return TRUE;
	}
	index--;
	return FALSE;
    }

    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    if (find_realized_shape_summary_in_node(group->getChild(i),
						    index, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLDatabaseSource::getRealizedShapeSummaryCount(void) const
{
    int count = this->d->compactIndexActive && this->d->compactIndex ?
	static_cast<int>(this->d->compactIndex->entries.size()) : 0;
    for (int i = 0; i < this->getNumChildren(); i++) {
	count += count_shapes_in_node(this->getChild(i));
	count += count_meshes_in_node(this->getChild(i));
    }
    return count;
}

SbBool
SoBRLDatabaseSource::getRealizedShapeSummary(int index,
	BObolRealizedShapeSummary &summary) const
{
    summary = BObolRealizedShapeSummary();
    if (index < 0)
	return FALSE;

    if (this->d->compactIndexActive && this->d->compactIndex) {
	const size_t compactCount = this->d->compactIndex->entries.size();
	if (static_cast<size_t>(index) < compactCount) {
	    summary = this->d->compactIndex->entries[
		static_cast<size_t>(index)].shapeSummary;
	    realized_shape_summary_owner(this, -1, summary);
	    return summary.valid;
	}
	index -= static_cast<int>(compactCount);
    }

    int remaining = index;
    for (int i = 0; i < this->getNumChildren(); i++) {
	if (find_realized_shape_summary_in_node(this->getChild(i),
						remaining, summary)) {
	    realized_shape_summary_owner(this, -1, summary);
	    return TRUE;
	}
    }

    return FALSE;
}

SbBool
SoBRLDatabaseSource::getSummary(BObolDatabaseSourceSummary &summary) const
{
    summary = BObolDatabaseSourceSummary();
    summary.valid = TRUE;
    summary.path = this->path.getValue();
    summary.instanceKey = source_effective_instance_key(this);
    summary.parentInstanceKey = this->parentInstanceKey.getValue();
    summary.occurrenceIndex = this->occurrenceIndex.getValue();
    summary.booleanOperation = this->booleanOperation.getValue();
    summary.displayName = this->displayName.getValue();
    summary.representationKey =
	source_effective_representation_key(this);
    summary.representationMode = this->representationMode.getValue();
    summary.drawMode = this->drawMode.getValue();
    summary.sourceRevision = this->sourceRevision.getValue();
    summary.inputsRevision = this->inputsRevision.getValue();
    summary.viewRevision = this->viewRevision.getValue();
    summary.realizedRevision = this->realizedRevision.getValue();
    summary.realizedSourceRevision = this->realizedSourceRevision.getValue();
    summary.realizedInputsRevision = this->realizedInputsRevision.getValue();
    summary.realizedViewRevision = this->realizedViewRevision.getValue();
    summary.realizationStatus = this->realizationStatus.getValue();
    summary.realizationDiagnostic = this->realizationDiagnostic.getValue();
    summary.realizationIdentity = this->realizationIdentity.getValue();
    summary.realizationRoleFlags = this->realizationRoleFlags.getValue();
    summary.realizationViewDependent =
	this->realizationViewDependent.getValue();
    summary.realizationCsgLodEnabled =
	this->realizationCsgLodEnabled.getValue();
    summary.realizationMeshLodEnabled =
	this->realizationMeshLodEnabled.getValue();
    summary.realizationViewScale = this->realizationViewScale.getValue();
    summary.realizationLodScale = this->realizationLodScale.getValue();
    summary.realizationViewWidth = this->realizationViewWidth.getValue();
    summary.realizationViewHeight = this->realizationViewHeight.getValue();
    summary.realizationBotThreshold =
	this->realizationBotThreshold.getValue();
    summary.realizationCurveScale = this->realizationCurveScale.getValue();
    summary.realizationPointScale = this->realizationPointScale.getValue();
    summary.visible = this->visible.getValue();
    summary.selected = this->selected.getValue();
    summary.highlighted = this->highlighted.getValue();
    summary.lineStyle = this->lineStyle.getValue();
    summary.lineWidth = this->lineWidth.getValue();
    summary.transparency = this->transparency.getValue();
    summary.materialColorValid = this->materialColorValid.getValue();
    summary.materialColor = this->materialColor.getValue();
    summary.materialRevision = this->materialRevision.getValue();
    summary.materialPolicy = this->materialPolicy.getValue();
    summary.databaseMetadataValid = this->databaseMetadataValid.getValue();
    summary.databaseRegionId = this->databaseRegionId.getValue();
    summary.databaseAirCode = this->databaseAirCode.getValue();
    summary.databaseMaterialId = this->databaseMaterialId.getValue();
    summary.databaseLos = this->databaseLos.getValue();
    summary.databaseMaterialColorValid =
	this->databaseMaterialColorValid.getValue();
    summary.databaseMaterialColor = this->databaseMaterialColor.getValue();
    summary.databaseMaterialShader = this->databaseMaterialShader.getValue();
    summary.colorOverride = this->colorOverride.getValue();
    summary.color = this->color.getValue();
    summary.selectedColor = this->selectedColor.getValue();
    summary.highlightedColor = this->highlightedColor.getValue();
    summary.ghostedColor = this->ghostedColor.getValue();
    summary.drawMatrixValid = this->drawMatrixValid.getValue();
    summary.drawMatrix = this->drawMatrix.getValue();
    summary.drawCenterValid = this->drawCenterValid.getValue();
    summary.drawCenter = this->drawCenter.getValue();
    summary.drawSizeValid = this->drawSizeValid.getValue();
    summary.drawSize = this->drawSize.getValue();
    summary.sourceBoundsValid = this->getEffectiveSourceBounds(
				    summary.sourceBounds);
    summary.stale = this->stale.getValue();
    summary.staleReason = this->staleReason.getValue();
    summary.realizedShapeCount = this->getRealizedShapeCount();
    summary.realizedMeshCount = this->getRealizedMeshCount();
    summary.realizedMaterialObjectCount =
	this->getRealizedMaterialObjectCount();
    return TRUE;
}
