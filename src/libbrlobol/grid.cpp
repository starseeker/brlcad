/*                         G R I D . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "bv.h"
#include "brlobol/grid.h"
#include "brlobol/lod_realization.h"
#include "brlobol/vlist_shape.h"

#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>
#include <Inventor/actions/SoSearchAction.h>
#include <Inventor/SoPath.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

SO_NODE_SOURCE(SoBRLGrid);

SoBRLGrid::SoBRLGrid(void)
{
    SO_NODE_CONSTRUCTOR(SoBRLGrid);

    SO_NODE_ADD_FIELD(overlayId, ("overlay::grid"));
    SO_NODE_ADD_FIELD(center, (0.0f, 0.0f, 0.0f));
    SO_NODE_ADD_FIELD(spacing, (1.0f));
    SO_NODE_ADD_FIELD(spacingV, (1.0f));
    SO_NODE_ADD_FIELD(divisions, (4));
    SO_NODE_ADD_FIELD(majorDivisionsH, (5));
    SO_NODE_ADD_FIELD(majorDivisionsV, (5));
    SO_NODE_ADD_FIELD(adaptive, (TRUE));
    SO_NODE_ADD_FIELD(visible, (TRUE));
    SO_NODE_ADD_FIELD(snapEnabled, (FALSE));
    SO_NODE_ADD_FIELD(overlayIntent, (TRUE));
    SO_NODE_ADD_FIELD(modelToView, (SbMatrix::identity()));
    SO_NODE_ADD_FIELD(viewScale, (1.0f));
    SO_NODE_ADD_FIELD(baseToLocal, (1.0f));
    SO_NODE_ADD_FIELD(aspect, (1.0f));
    SO_NODE_ADD_FIELD(viewportWidth, (512));
    SO_NODE_ADD_FIELD(viewportHeight, (512));
    SO_NODE_ADD_FIELD(targetPixelSpacing, (56.0f));
    SO_NODE_ADD_FIELD(minimumPixelSpacing, (8.0f));
    SO_NODE_ADD_FIELD(minorColor, (SbColor(0.24f, 0.24f, 0.24f)));
    SO_NODE_ADD_FIELD(majorColor, (SbColor(0.42f, 0.42f, 0.42f)));
    SO_NODE_ADD_FIELD(axisColor, (SbColor(0.58f, 0.58f, 0.58f)));
    SO_NODE_ADD_FIELD(effectiveSpacingH, (1.0f));
    SO_NODE_ADD_FIELD(effectiveSpacingV, (1.0f));
    SO_NODE_ADD_FIELD(pixelSpacingH, (1.0f));
    SO_NODE_ADD_FIELD(pixelSpacingV, (1.0f));
    SO_NODE_ADD_FIELD(minorSegmentCount, (0));
    SO_NODE_ADD_FIELD(majorSegmentCount, (0));
    SO_NODE_ADD_FIELD(axisSegmentCount, (0));
}

SoBRLGrid::~SoBRLGrid(void)
{
}

void
SoBRLGrid::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLGrid, SoSeparator, "Separator");
}

static SbMatrix
grid_sbmatrix_from_brl_mat(const mat_t mat)
{
    if (!mat)
	return SbMatrix::identity();

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

static SbColor
grid_color_from_int_rgb(const int rgb[3], float scale)
{
    if (!rgb)
	return SbColor(scale, scale, scale);

    const float r = std::max(0.0f, std::min(1.0f,
					    static_cast<float>(rgb[0]) / 255.0f * scale));
    const float g = std::max(0.0f, std::min(1.0f,
					    static_cast<float>(rgb[1]) / 255.0f * scale));
    const float b = std::max(0.0f, std::min(1.0f,
					    static_cast<float>(rgb[2]) / 255.0f * scale));
    return SbColor(r, g, b);
}

static double
grid_nice_spacing(double configuredSpacing, double pixelSize,
		  double targetPixels, double minimumPixels, SbBool adaptive)
{
    if (configuredSpacing <= 0.0)
	configuredSpacing = 1.0;
    if (pixelSize <= 0.0)
	return configuredSpacing;

    const double configuredPixels = configuredSpacing / pixelSize;
    if (!adaptive && configuredPixels >= minimumPixels)
	return configuredSpacing;

    const double desiredSpacing = pixelSize * std::max(targetPixels,
				  minimumPixels);
    if (desiredSpacing <= configuredSpacing)
	return configuredSpacing;

    const double ratio = desiredSpacing / configuredSpacing;
    const double decade = std::pow(10.0, std::floor(std::log10(ratio)));
    const double candidates[4] = {1.0, 2.0, 5.0, 10.0};
    for (size_t i = 0; i < 4; i++) {
	const double factor = candidates[i] * decade;
	if (ratio <= factor)
	    return configuredSpacing * factor;
    }

    return configuredSpacing * 10.0 * decade;
}

static void
grid_append_line(std::vector<SbVec3f> &points,
		 std::vector<int32_t> &commands,
		 float x1,
		 float y1,
		 float x2,
		 float y2)
{
    points.push_back(SbVec3f(x1, y1, 0.0f));
    commands.push_back(SoBRLVListShape::MOVE);
    points.push_back(SbVec3f(x2, y2, 0.0f));
    commands.push_back(SoBRLVListShape::DRAW);
}

static SoBRLVListShape *
grid_make_shape(const SbString &overlayId,
		const char *suffix,
		const char *geometryName,
		const SbColor &color,
		int lineWidth,
		const std::vector<SbVec3f> &points,
		const std::vector<int32_t> &commands)
{
    if (points.empty() || commands.empty())
	return NULL;

    SbString path = overlayId;
    path += suffix ? suffix : "";

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->sourcePath = path;
    shape->displayName = overlayId;
    shape->geometryName = geometryName ? geometryName : "grid";
    shape->sourceIdentity = path;
    shape->cacheIdentity = path;
    shape->databaseIntent = FALSE;
    shape->overlayIntent = TRUE;
    shape->hudIntent = TRUE;
    shape->localSource = TRUE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = BRLOBOL_LOD_DRAW_DIAGNOSTIC;
    shape->recordRole = "overlay";
    shape->geometryKind = "line";
    shape->visible = TRUE;
    shape->selectable = FALSE;
    shape->colorOverride = TRUE;
    shape->color = color;
    shape->lineWidth = lineWidth > 0 ? lineWidth : 1;
    shape->setLineSet(points.data(), commands.data(),
		      static_cast<int>(points.size()));
    return shape;
}

static int
grid_line_class(int index, int majorDivisions)
{
    if (index == 0)
	return 2;
    if (majorDivisions > 0 && index % majorDivisions == 0)
	return 1;
    return 0;
}

static void
grid_anchor_view_local(const SoBRLGrid *grid, double sf, point_t out)
{
    const SbVec3f c = grid->center.getValue();
    SbVec3f viewAnchor;
    grid->modelToView.getValue().multVecMatrix(c, viewAnchor);
    VSET(out, viewAnchor[0] * sf, viewAnchor[1] * sf, viewAnchor[2] * sf);
}

SoBRLVListShape *
SoBRLGrid::rebuildGeometry(void)
{
    this->removeAllChildren();
    this->minorSegmentCount = 0;
    this->majorSegmentCount = 0;
    this->axisSegmentCount = 0;
    if (!this->visible.getValue())
	return NULL;

    int width = this->viewportWidth.getValue();
    int height = this->viewportHeight.getValue();
    if (width <= 0)
	width = 512;
    if (height <= 0)
	height = 512;

    double aspectValue = this->aspect.getValue();
    if (aspectValue <= 0.0)
	aspectValue = static_cast<double>(width) / static_cast<double>(height);
    if (aspectValue <= 0.0)
	aspectValue = 1.0;

    double sf = static_cast<double>(this->viewScale.getValue()) *
		static_cast<double>(this->baseToLocal.getValue());
    if (sf <= 0.0)
	sf = 1.0;

    const double pixelSize = 2.0 * sf / static_cast<double>(width);
    const double configuredH = std::max(1.0e-12,
					static_cast<double>(this->spacing.getValue()) *
					static_cast<double>(this->baseToLocal.getValue()));
    const double configuredV = std::max(1.0e-12,
					static_cast<double>(this->spacingV.getValue()) *
					static_cast<double>(this->baseToLocal.getValue()));
    const double targetPixels = std::max(1.0,
					 static_cast<double>(this->targetPixelSpacing.getValue()));
    const double minimumPixels = std::max(1.0,
					  static_cast<double>(this->minimumPixelSpacing.getValue()));
    const double effectiveH = grid_nice_spacing(configuredH, pixelSize,
			      targetPixels, minimumPixels,
			      this->adaptive.getValue());
    const double effectiveV = grid_nice_spacing(configuredV, pixelSize,
			      targetPixels, minimumPixels,
			      this->adaptive.getValue());
    this->effectiveSpacingH = static_cast<float>(effectiveH);
    this->effectiveSpacingV = static_cast<float>(effectiveV);
    this->pixelSpacingH = static_cast<float>(effectiveH / pixelSize);
    this->pixelSpacingV = static_cast<float>(effectiveV / pixelSize);

    point_t anchorLocal;
    grid_anchor_view_local(this, sf, anchorLocal);

    const double left = -sf;
    const double right = sf;
    const double bottom = -sf / aspectValue;
    const double top = sf / aspectValue;
    const int majorH = std::max(1, this->majorDivisionsH.getValue());
    const int majorV = std::max(1, this->majorDivisionsV.getValue());

    std::vector<SbVec3f> minorPoints;
    std::vector<int32_t> minorCommands;
    std::vector<SbVec3f> majorPoints;
    std::vector<int32_t> majorCommands;
    std::vector<SbVec3f> axisPoints;
    std::vector<int32_t> axisCommands;

    const int firstX = static_cast<int>(std::floor((left - anchorLocal[X]) /
					effectiveH));
    const int lastX = static_cast<int>(std::ceil((right - anchorLocal[X]) /
				       effectiveH));
    const int firstY = static_cast<int>(std::floor((bottom - anchorLocal[Y]) /
					effectiveV));
    const int lastY = static_cast<int>(std::ceil((top - anchorLocal[Y]) /
				       effectiveV));

    for (int i = firstX; i <= lastX; i++) {
	const double localX = anchorLocal[X] + static_cast<double>(i) *
			      effectiveH;
	const float px = static_cast<float>((localX / sf + 1.0) * 0.5 *
					    static_cast<double>(width));
	if (px < -1.0f || px > static_cast<float>(width) + 1.0f)
	    continue;
	const int klass = grid_line_class(i, majorH);
	if (klass == 2)
	    grid_append_line(axisPoints, axisCommands, px, 0.0f, px,
			     static_cast<float>(height));
	else if (klass == 1)
	    grid_append_line(majorPoints, majorCommands, px, 0.0f, px,
			     static_cast<float>(height));
	else if (this->pixelSpacingH.getValue() >= minimumPixels)
	    grid_append_line(minorPoints, minorCommands, px, 0.0f, px,
			     static_cast<float>(height));
    }

    for (int i = firstY; i <= lastY; i++) {
	const double localY = anchorLocal[Y] + static_cast<double>(i) *
			      effectiveV;
	const float py = static_cast<float>(((localY / sf) * aspectValue +
					     1.0) * 0.5 *
					    static_cast<double>(height));
	if (py < -1.0f || py > static_cast<float>(height) + 1.0f)
	    continue;
	const int klass = grid_line_class(i, majorV);
	if (klass == 2)
	    grid_append_line(axisPoints, axisCommands, 0.0f, py,
			     static_cast<float>(width), py);
	else if (klass == 1)
	    grid_append_line(majorPoints, majorCommands, 0.0f, py,
			     static_cast<float>(width), py);
	else if (this->pixelSpacingV.getValue() >= minimumPixels)
	    grid_append_line(minorPoints, minorCommands, 0.0f, py,
			     static_cast<float>(width), py);
    }

    SoBRLVListShape *minor = grid_make_shape(this->overlayId.getValue(),
			     "::minor", "grid-minor", this->minorColor.getValue(),
			     1, minorPoints, minorCommands);
    SoBRLVListShape *major = grid_make_shape(this->overlayId.getValue(),
			     "::major", "grid-major", this->majorColor.getValue(),
			     1, majorPoints, majorCommands);
    SoBRLVListShape *axis = grid_make_shape(this->overlayId.getValue(),
					    "::axis", "grid-axis", this->axisColor.getValue(),
					    2, axisPoints, axisCommands);
    this->minorSegmentCount = minor ? minor->getSegmentCount() : 0;
    this->majorSegmentCount = major ? major->getSegmentCount() : 0;
    this->axisSegmentCount = axis ? axis->getSegmentCount() : 0;
    if (!minor && !major && !axis)
	return NULL;

    SoHUDKit *hud = new SoHUDKit;
    if (minor)
	hud->addWidget(minor);
    if (major)
	hud->addWidget(major);
    if (axis)
	hud->addWidget(axis);
    this->addChild(hud);

    return axis ? axis : (major ? major : minor);
}

static SoBRLVListShape *
grid_find_shape(const SoBRLGrid *grid, const char *geometryName)
{
    if (!grid || !geometryName)
	return NULL;

    SoSearchAction search;
    search.setType(SoBRLVListShape::getClassTypeId());
    search.setInterest(SoSearchAction::ALL);
    search.setSearchingAll(TRUE);
    search.apply(const_cast<SoBRLGrid *>(grid));
    SoPathList &paths = search.getPaths();
    for (int i = 0; i < paths.getLength(); i++) {
	SoPath *path = paths[i];
	if (!path)
	    continue;
	SoNode *tail = path->getTail();
	if (!tail || !tail->isOfType(SoBRLVListShape::getClassTypeId()))
	    continue;
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(tail);
	if (bu_strcmp(shape->geometryName.getValue().getString(),
		   geometryName) == 0)
	    return shape;
    }
    return NULL;
}

SoBRLVListShape *
SoBRLGrid::getGeometryShape(void) const
{
    SoBRLVListShape *axis = this->getAxisShape();
    if (axis)
	return axis;
    SoBRLVListShape *major = this->getMajorShape();
    return major ? major : this->getMinorShape();
}

SoBRLVListShape *
SoBRLGrid::getMinorShape(void) const
{
    return grid_find_shape(this, "grid-minor");
}

SoBRLVListShape *
SoBRLGrid::getMajorShape(void) const
{
    return grid_find_shape(this, "grid-major");
}

SoBRLVListShape *
SoBRLGrid::getAxisShape(void) const
{
    return grid_find_shape(this, "grid-axis");
}

int
SoBRLGrid::getMinorSegmentCount(void) const
{
    SoBRLVListShape *shape = this->getMinorShape();
    return shape ? shape->getSegmentCount() : this->minorSegmentCount.getValue();
}

int
SoBRLGrid::getMajorSegmentCount(void) const
{
    SoBRLVListShape *shape = this->getMajorShape();
    return shape ? shape->getSegmentCount() : this->majorSegmentCount.getValue();
}

int
SoBRLGrid::getAxisSegmentCount(void) const
{
    SoBRLVListShape *shape = this->getAxisShape();
    return shape ? shape->getSegmentCount() : this->axisSegmentCount.getValue();
}

int
SoBRLGrid::getTotalSegmentCount(void) const
{
    return this->getMinorSegmentCount() + this->getMajorSegmentCount() +
	   this->getAxisSegmentCount();
}

int
brlobol_grid_configure_from_view(SoBRLGrid *grid,
				 const struct bv_grid_state *state,
				 const mat_t model2view,
				 fastf_t view_scale,
				 fastf_t base2local,
				 int width,
				 int height)
{
    if (!grid || !state)
	return 0;

    grid->visible = state->draw ? TRUE : FALSE;
    grid->snapEnabled = state->snap ? TRUE : FALSE;
    grid->overlayIntent = TRUE;
    grid->center = SbVec3f(static_cast<float>(state->anchor[X]),
			   static_cast<float>(state->anchor[Y]),
			   static_cast<float>(state->anchor[Z]));
    grid->spacing = static_cast<float>(state->res_h > SMALL_FASTF ?
				       state->res_h : 1.0);
    grid->spacingV = static_cast<float>(state->res_v > SMALL_FASTF ?
					state->res_v : grid->spacing.getValue());
    grid->divisions = std::max(state->res_major_h, state->res_major_v);
    grid->majorDivisionsH = state->res_major_h > 0 ? state->res_major_h : 5;
    grid->majorDivisionsV = state->res_major_v > 0 ? state->res_major_v : 5;
    grid->adaptive = state->adaptive ? TRUE : FALSE;
    grid->modelToView = grid_sbmatrix_from_brl_mat(model2view);
    grid->viewScale = static_cast<float>(view_scale > SMALL_FASTF ?
					 view_scale : 1.0);
    grid->baseToLocal = static_cast<float>(base2local > SMALL_FASTF ?
					   base2local : 1.0);
    grid->viewportWidth = width > 0 ? width : 512;
    grid->viewportHeight = height > 0 ? height : 512;
    grid->aspect = width > 0 && height > 0 ?
		   static_cast<float>(width) / static_cast<float>(height) : 1.0f;
    grid->minorColor = grid_color_from_int_rgb(state->color, 0.24f);
    grid->majorColor = grid_color_from_int_rgb(state->color, 0.42f);
    grid->axisColor = grid_color_from_int_rgb(state->color, 0.62f);
    grid->rebuildGeometry();
    return (state->draw || state->snap) ? 1 : 0;
}

int
brlobol_grid_configure_from_view_context(SoBRLGrid *grid,
	const struct bv_grid_state *state,
	const void *view_ctx)
{
    const struct bv *view =
	bv_context_view_const((const struct bv_context *)view_ctx);
    if (!grid || !state || !view)
	return 0;

    mat_t model2view;
    MAT_IDN(model2view);
    (void)bv_model2view_get(model2view, view);
    return brlobol_grid_configure_from_view(grid, state, model2view,
					    bv_scale_get(view),
					    bv_base2local_get(view),
					    bv_width_get(view),
					    bv_height_get(view));
}
