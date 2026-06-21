/*                E D I T _ P R E V I E W . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/edit_preview.h"
#include "brlobol/vlist_shape.h"

#include <Inventor/nodes/SoMatrixTransform.h>
#include <Inventor/sensors/SoFieldSensor.h>

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

void
SoBRLEditPreview::clearPreview(void)
{
    this->removeAllChildren();
    this->realizedSourceRevision = this->sourceRevision.getValue();
    this->realizedInputsRevision = this->inputsRevision.getValue();
    this->stale = FALSE;
    this->previewStatus = EMPTY;
}

SoBRLVListShape *
SoBRLEditPreview::appendLineSet(const SbString &identity,
	const SbVec3f *points,
	const int32_t *commands,
	int count)
{
    if (!points || !commands || count <= 0) {
	this->previewStatus = FAILED;
	this->stale = TRUE;
	return NULL;
    }

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->sourcePath = identity.getLength() ? identity : this->previewId.getValue();
    shape->sourceName = this->previewId.getValue();
    shape->sourceType = "edit-preview";
    shape->sourceId = this->sourceRevision.getValue();
    shape->editEmphasis = TRUE;
    const SbString &intentId = this->editIntentId.getValue();
    const SbString &intentRole = this->editIntentRole.getValue();
    shape->editIntentId = intentId.getLength() ? intentId : this->previewId.getValue();
    shape->editIntentRole = intentRole.getLength() ? intentRole : SbString("preview");
    shape->setLineSet(points, commands, count);
    return shape;
}

void
SoBRLEditPreview::markCurrent(void)
{
    this->realizedSourceRevision = this->sourceRevision.getValue();
    this->realizedInputsRevision = this->inputsRevision.getValue();
    this->stale = FALSE;
    this->previewStatus = CURRENT;
}

SoBRLVListShape *
SoBRLEditPreview::setLineSet(const SbString &identity,
	const SbVec3f *points,
	const int32_t *commands,
	int count)
{
    this->removeAllChildren();

    SoBRLVListShape *shape = this->appendLineSet(identity, points, commands, count);
    if (!shape)
	return NULL;

    this->addChild(shape);
    this->markCurrent();
    return shape;
}

SoBRLVListShape *
SoBRLEditPreview::setTransformedLineSet(const SbString &identity,
	const SbMatrix &matrix,
	const SbVec3f *points,
	const int32_t *commands,
	int count)
{
    this->removeAllChildren();

    SoBRLVListShape *shape = this->appendLineSet(identity, points, commands, count);
    if (!shape)
	return NULL;

    SoSeparator *root = new SoSeparator;
    SoMatrixTransform *transform = new SoMatrixTransform;
    transform->matrix = matrix;
    root->addChild(transform);
    root->addChild(shape);
    this->addChild(root);
    this->markCurrent();
    return shape;
}
