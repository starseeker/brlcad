/*                E D I T _ P R E V I E W . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BEditPreview.h"
#include "BObol/BLodRealization.h"
#include "BObol/BVListShape.h"
#include "scalar_publication_private.h"

#include <Inventor/misc/SoChildList.h>
#include <Inventor/nodes/SoMatrixTransform.h>
#include <Inventor/sensors/SoFieldSensor.h>
#include <Inventor/tools/SbModernUtils.h>

#include <exception>
#include <vector>

SO_NODE_SOURCE(SoBRLEditPreview);

SoBRLEditPreview::SoBRLEditPreview(void) :
    previewIdSensor(NULL),
    editIntentIdSensor(NULL),
    editIntentRoleSensor(NULL),
    sourceRevisionSensor(NULL),
    inputsRevisionSensor(NULL)
{
    SO_NODE_CONSTRUCTOR(SoBRLEditPreview);

    SO_NODE_DEFINE_ENUM_VALUE(PreviewStatus, EMPTY);
    SO_NODE_DEFINE_ENUM_VALUE(PreviewStatus, CURRENT);
    SO_NODE_DEFINE_ENUM_VALUE(PreviewStatus, STALE);
    SO_NODE_DEFINE_ENUM_VALUE(PreviewStatus, FAILED);

    SO_NODE_ADD_FIELD(previewId, (""));
    SO_NODE_ADD_FIELD(editIntentId, (""));
    SO_NODE_ADD_FIELD(editIntentRole, ("preview"));
    SO_NODE_ADD_FIELD(sourceRevision, (0));
    SO_NODE_ADD_FIELD(inputsRevision, (0));
    SO_NODE_ADD_FIELD(realizedSourceRevision, (0));
    SO_NODE_ADD_FIELD(realizedInputsRevision, (0));
    SO_NODE_ADD_FIELD(previewStatus, (EMPTY));
    SO_NODE_SET_SF_ENUM_TYPE(previewStatus, PreviewStatus);
    SO_NODE_ADD_FIELD(stale, (FALSE));

    this->attachFieldSensors();
}

SoBRLEditPreview::~SoBRLEditPreview(void)
{
    this->detachFieldSensors();
}

void
SoBRLEditPreview::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLEditPreview, SoSeparator, "Separator");
}

void
SoBRLEditPreview::fieldSensorCB(void *data, SoSensor *UNUSED(sensor))
{
    SoBRLEditPreview *preview = static_cast<SoBRLEditPreview *>(data);
    if (preview)
	preview->markStale();
}

void
SoBRLEditPreview::attachFieldSensors(void)
{
    this->previewIdSensor = new SoFieldSensor(SoBRLEditPreview::fieldSensorCB, this);
    this->previewIdSensor->setPriority(0);
    this->previewIdSensor->attach(&this->previewId);

    this->editIntentIdSensor = new SoFieldSensor(SoBRLEditPreview::fieldSensorCB, this);
    this->editIntentIdSensor->setPriority(0);
    this->editIntentIdSensor->attach(&this->editIntentId);

    this->editIntentRoleSensor = new SoFieldSensor(SoBRLEditPreview::fieldSensorCB, this);
    this->editIntentRoleSensor->setPriority(0);
    this->editIntentRoleSensor->attach(&this->editIntentRole);

    this->sourceRevisionSensor = new SoFieldSensor(SoBRLEditPreview::fieldSensorCB, this);
    this->sourceRevisionSensor->setPriority(0);
    this->sourceRevisionSensor->attach(&this->sourceRevision);

    this->inputsRevisionSensor = new SoFieldSensor(SoBRLEditPreview::fieldSensorCB, this);
    this->inputsRevisionSensor->setPriority(0);
    this->inputsRevisionSensor->attach(&this->inputsRevision);
}

void
SoBRLEditPreview::detachFieldSensors(void)
{
    delete this->previewIdSensor;
    this->previewIdSensor = NULL;
    delete this->editIntentIdSensor;
    this->editIntentIdSensor = NULL;
    delete this->editIntentRoleSensor;
    this->editIntentRoleSensor = NULL;
    delete this->sourceRevisionSensor;
    this->sourceRevisionSensor = NULL;
    delete this->inputsRevisionSensor;
    this->inputsRevisionSensor = NULL;
}

void
SoBRLEditPreview::markStale(void)
{
    this->stale = TRUE;
    this->previewStatus = STALE;
}

void
SoBRLEditPreview::setEditIntent(const SbString &id, const SbString &role)
{
    this->editIntentId = id;
    this->editIntentRole = role.getLength() ? role : SbString("preview");
    this->markStale();
}

void
SoBRLEditPreview::markSourceRevision(uint32_t revision)
{
    this->sourceRevision = revision;
    this->markStale();
}

void
SoBRLEditPreview::markInputsRevision(uint32_t revision)
{
    this->inputsRevision = revision;
    this->markStale();
}

SbBool
SoBRLEditPreview::needsRealization(void) const
{
    return this->stale.getValue() ||
	   this->realizedSourceRevision.getValue() != this->sourceRevision.getValue() ||
	   this->realizedInputsRevision.getValue() != this->inputsRevision.getValue();
}

static void
edit_preview_publish(SoBRLEditPreview &preview, SoNode *child,
	uint32_t realizedSourceRevision, uint32_t realizedInputsRevision,
	SoBRLEditPreview::PreviewStatus status, SbBool stale)
{
    std::vector<SoNode *> publishedChildren;
    if (child)
	publishedChildren.push_back(child);
    auto replacement =
	preview.getChildren()->prepareReplacement(publishedChildren);
    PreparedFieldNotifications<4> notifications(preview, {{
	{&preview.realizedSourceRevision,
	    preview.realizedSourceRevision.getValue() != realizedSourceRevision},
	{&preview.realizedInputsRevision,
	    preview.realizedInputsRevision.getValue() != realizedInputsRevision},
	{&preview.previewStatus,
	    preview.previewStatus.getValue() != static_cast<int>(status)},
	{&preview.stale, preview.stale.getValue() != stale}
    }});

    replacement->commit();
    preview.realizedSourceRevision = realizedSourceRevision;
    preview.realizedInputsRevision = realizedInputsRevision;
    preview.previewStatus = status;
    preview.stale = stale;

    notifications.restore();
    std::exception_ptr failure;
    try { replacement->notify(); }
    catch (...) { failure = std::current_exception(); }
    notifications.notify(failure);
    if (failure)
	std::rethrow_exception(failure);
}

static SbModernUtils::SoNodeRef
edit_preview_make_line_shape(const SoBRLEditPreview &preview,
	const SbString &identity, const SbVec3f *points,
	const int32_t *commands, int count)
{
    if (!points || !commands || count <= 0)
	return SbModernUtils::SoNodeRef(nullptr);

    SbModernUtils::SoNodeRef owner(new SoBRLVListShape);
    auto *shape = static_cast<SoBRLVListShape *>(owner.get());
    shape->sourcePath = identity.getLength() ? identity : preview.previewId.getValue();
    shape->sourceName = preview.previewId.getValue();
    shape->sourceType = "edit-preview";
    shape->sourceId = preview.sourceRevision.getValue();
    shape->displayName = preview.previewId.getValue();
    shape->geometryName = preview.previewId.getValue();
    shape->sourceIdentity = shape->sourcePath.getValue();
    shape->cacheIdentity = shape->sourcePath.getValue();
    shape->databaseIntent = FALSE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = TRUE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = BOBOL_LOD_DRAW_DIAGNOSTIC;
    shape->recordRole = "edit-preview";
    shape->geometryKind = "line";
    shape->editEmphasis = TRUE;
    /* Edit previews must render in the edit-emphasis color, not the source
     * solid's material.  Without an explicit override the shape inherits the
     * edited solid's material (e.g. LIGHT's red) via cad_shape_color().  Apply
     * the conventional edit color (yellow, matching MGED's cs_edit_info
     * default); a caller-supplied style color (store_apply_vlist_style) still
     * overrides this if the edit transaction ever plumbs a scheme color. */
    shape->colorOverride = TRUE;
    shape->color = SbColor(1.0f, 1.0f, 0.0f);
    const SbString &intentId = preview.editIntentId.getValue();
    const SbString &intentRole = preview.editIntentRole.getValue();
    shape->editIntentId = intentId.getLength() ? intentId : preview.previewId.getValue();
    shape->editIntentRole = intentRole.getLength() ? intentRole : SbString("preview");
    shape->setLineSet(points, commands, count);
    return owner;
}

static SoBRLVListShape *
edit_preview_publish_failure(SoBRLEditPreview &preview)
{
    if (!preview.getNumChildren() && preview.stale.getValue() &&
	preview.previewStatus.getValue() == SoBRLEditPreview::FAILED)
	return NULL;
    edit_preview_publish(preview, NULL,
	preview.realizedSourceRevision.getValue(),
	preview.realizedInputsRevision.getValue(), SoBRLEditPreview::FAILED, TRUE);
    return NULL;
}

void
SoBRLEditPreview::clearPreview(void)
{
    const uint32_t source = this->sourceRevision.getValue();
    const uint32_t inputs = this->inputsRevision.getValue();
    if (!this->getNumChildren() &&
	this->realizedSourceRevision.getValue() == source &&
	this->realizedInputsRevision.getValue() == inputs &&
	!this->stale.getValue() &&
	this->previewStatus.getValue() == EMPTY)
	return;
    edit_preview_publish(*this, NULL, source, inputs, EMPTY, FALSE);
}

SoBRLVListShape *
SoBRLEditPreview::setLineSet(const SbString &identity,
			     const SbVec3f *points,
			     const int32_t *commands,
			     int count)
{
    SbModernUtils::SoNodeRef owner = edit_preview_make_line_shape(
	*this, identity, points, commands, count);
    if (!owner)
	return edit_preview_publish_failure(*this);
    auto *shape = static_cast<SoBRLVListShape *>(owner.get());
    edit_preview_publish(*this, shape, this->sourceRevision.getValue(),
	this->inputsRevision.getValue(), CURRENT, FALSE);
    return shape;
}

SoBRLVListShape *
SoBRLEditPreview::setTransformedLineSet(const SbString &identity,
					const SbMatrix &matrix,
					const SbVec3f *points,
					const int32_t *commands,
					int count)
{
    SbModernUtils::SoNodeRef shapeOwner = edit_preview_make_line_shape(
	*this, identity, points, commands, count);
    if (!shapeOwner)
	return edit_preview_publish_failure(*this);
    auto *shape = static_cast<SoBRLVListShape *>(shapeOwner.get());

    SbModernUtils::SoNodeRef rootOwner(new SoSeparator);
    auto *root = static_cast<SoSeparator *>(rootOwner.get());
    SbModernUtils::SoNodeRef transformOwner(new SoMatrixTransform);
    auto *transform = static_cast<SoMatrixTransform *>(transformOwner.get());
    transform->matrix = matrix;
    root->addChild(transform);
    root->addChild(shape);
    edit_preview_publish(*this, root, this->sourceRevision.getValue(),
	this->inputsRevision.getValue(), CURRENT, FALSE);
    return shape;
}
