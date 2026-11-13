/* BRL-CAD
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
#ifndef BOBOL_SCALAR_PUBLICATION_PRIVATE_H
#define BOBOL_SCALAR_PUBLICATION_PRIVATE_H

#include <Inventor/fields/SoSField.h>
#include <Inventor/fields/SoSFNode.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFInt32.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFVec2f.h>
#include <Inventor/fields/SoSFVec3f.h>
#include <Inventor/fields/SoSFColor.h>
#include <Inventor/fields/SoSFMatrix.h>
#include <Inventor/fields/SoSFPlane.h>
#include <Inventor/fields/SoSFRotation.h>
#include <Inventor/fields/SoSFEnum.h>
#include <Inventor/fields/SoFieldData.h>
#include <Inventor/misc/SoChildList.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/sensors/SoFieldSensor.h>
#include <algorithm>
#include <array>
#include <exception>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

struct PublicationFieldChange {
    SoField *field;
    bool changed;
    SbBool notifications = FALSE;
    SoFieldSensor *handledSensor = nullptr;
};

/* A prepared field notification has already committed every participant in
 * its publication.  Source nodes consume this one-shot marker before
 * propagating the notification, which distinguishes that completed event
 * from a reentrant direct field write without suppressing the latter. */
class PreparedFieldNotificationScope {
public:
    explicit PreparedFieldNotificationScope(SoField &field) :
	previous(active())
    {
	active() = &field;
    }
    ~PreparedFieldNotificationScope() { active() = this->previous; }

    PreparedFieldNotificationScope(const PreparedFieldNotificationScope &) = delete;
    PreparedFieldNotificationScope &operator=(const PreparedFieldNotificationScope &) = delete;

    static bool consume(const SoField *field)
    {
	if (!field || active() != field)
	    return false;
	active() = nullptr;
	return true;
    }

private:
    static SoField *&active()
    {
	static thread_local SoField *field = nullptr;
	return field;
    }

    SoField *previous;
};

/* A publication may span several nodes. Restore every participant before
 * delivering callbacks, so an observer can safely initiate a later update. */
template <typename Fields>
class PreparedNotifications {
public:
    PreparedNotifications(SoNode &target, Fields changes) :
	node(target), fields(std::move(changes)), enabled(target.isNotifyEnabled())
    {
	for (auto &field : this->fields)
	    field.notifications = field.field->enableNotify(FALSE);
	this->node.enableNotify(FALSE);
    }
    ~PreparedNotifications() { this->restore(); }
    PreparedNotifications(const PreparedNotifications &) = delete;
    PreparedNotifications &operator=(const PreparedNotifications &) = delete;

    bool changed() const
    {
	return std::any_of(this->fields.begin(), this->fields.end(),
	    [](const PublicationFieldChange &field) { return field.changed; });
    }
    void restore()
    {
	if (this->restored)
	    return;
	for (const auto &field : this->fields)
	    field.field->enableNotify(field.notifications);
	this->node.enableNotify(this->enabled);
	this->restored = true;
    }
    void notify(std::exception_ptr &failure)
    {
	this->restore();
	bool fieldNotified = false;
	for (const auto &field : this->fields) {
	    if (!field.changed || !field.notifications)
		continue;
	    PreparedFieldNotificationScope prepared(*field.field);
	    try { field.field->touch(field.handledSensor); }
	    catch (...) { if (!failure) failure = std::current_exception(); }
	    fieldNotified = true;
	}
	if (!fieldNotified && this->enabled) {
	    try { this->node.touch(); }
	    catch (...) { if (!failure) failure = std::current_exception(); }
	}
    }

private:
    SoNode &node;
    Fields fields;
    SbBool enabled;
    bool restored = false;
};

template <size_t Count>
using PreparedFieldNotifications = PreparedNotifications<std::array<PublicationFieldChange, Count>>;

/* Only these scalar field types commit without allocating. Strings own a
 * prepared movable value. Node and multi-value fields require a dedicated
 * ownership transaction and are excluded here. */
static bool
publication_scalar_field(const SoField &field)
{
    const SoType type = field.getTypeId();
    if (!type.isDerivedFrom(SoSField::getClassTypeId()) || type == SoSFNode::getClassTypeId())
	return false;
    if (type == SoSFString::getClassTypeId() || type == SoSFInt32::getClassTypeId() ||
	type == SoSFUInt32::getClassTypeId() || type == SoSFBool::getClassTypeId() ||
	type == SoSFFloat::getClassTypeId() || type == SoSFVec2f::getClassTypeId() ||
	type == SoSFVec3f::getClassTypeId() ||
	type == SoSFColor::getClassTypeId() || type == SoSFMatrix::getClassTypeId() ||
	type == SoSFPlane::getClassTypeId() ||
	type == SoSFRotation::getClassTypeId() ||
	type.isDerivedFrom(SoSFEnum::getClassTypeId()))
	return true;
    throw std::logic_error("publication has unsupported scalar metadata");
}

class PreparedScalarValues {
public:
    void reserve(size_t count) { this->values.reserve(count); }
    bool prepare(SoField &field, const SoField *value)
    {
	if (!publication_scalar_field(field)) return false;
	if (!value || value->getTypeId() != field.getTypeId())
	    throw std::logic_error("publication has incompatible scalar metadata");
	if (field.isSame(*value)) return false;
	const bool string = field.getTypeId() == SoSFString::getClassTypeId();
	this->values.push_back({&field, value, string ?
	    static_cast<const SoSFString *>(value)->getValue() : SbString(), string,
	    field.getTypeId().isDerivedFrom(SoSFEnum::getClassTypeId()) != FALSE});
	return true;
    }
    void commit()
    {
	for (Scalar &scalar : this->values) {
	    if (scalar.string)
		static_cast<SoSFString *>(scalar.field)->setValue(std::move(scalar.text));
	    else if (scalar.enumeration)
		/* Enum assignment also copies mapping arrays. Only the prepared
		 * numeric value changes during this quiet commit. */
		static_cast<SoSFEnum *>(scalar.field)->setValue(
		    static_cast<const SoSFEnum *>(scalar.value)->getValue());
	    else
		scalar.field->copyFrom(*scalar.value);
	}
    }
private:
    struct Scalar {
	SoField *field;
	const SoField *value;
	SbString text;
	bool string;
	bool enumeration;
    };
    std::vector<Scalar> values;
};

inline void
copy_publication_scalar_fields(SoNode &target, const SoNode &source)
{
    const SoFieldData *data = static_cast<const SoFieldContainer &>(source).getFieldData();
    for (int i = 0; i < data->getNumFields(); ++i) {
	const SoField *field = data->getField(&source, i);
	if (!publication_scalar_field(*field)) continue;
	SoField *value = target.getField(data->getFieldName(i));
	if (!value || value->getTypeId() != field->getTypeId())
	    throw std::logic_error("publication has unsupported derived metadata");
	value->copyFrom(*field);
    }
}

class PreparedScalarFields {
public:
    PreparedScalarFields(SoNode &source, const SoNode &next,
	std::initializer_list<SoFieldSensor *> handledSensors = {})
    {
	const SoFieldData *data = static_cast<const SoFieldContainer &>(source).getFieldData();
	std::vector<PublicationFieldChange> changes;
	changes.reserve(size_t(data->getNumFields()));
	this->values.reserve(size_t(data->getNumFields()));
	for (int i = 0; i < data->getNumFields(); ++i) {
	    SoField *field = data->getField(&source, i);
	    SoFieldSensor *handled = nullptr;
	    for (SoFieldSensor *sensor : handledSensors)
		if (sensor && sensor->getAttachedField() == field) handled = sensor;
	    changes.push_back({field, this->values.prepare(*field, next.getField(data->getFieldName(i))),
		FALSE, handled});
	}
	this->notifications = std::make_unique<PreparedNotifications<std::vector<PublicationFieldChange>>>(
	    source, std::move(changes));
    }
    bool changed() const { return this->notifications->changed(); }
    void commit() { this->values.commit(); }
    void restore() { this->notifications->restore(); }
    void notify(std::exception_ptr &failure) { this->notifications->notify(failure); }
private:
    PreparedScalarValues values;
    std::unique_ptr<PreparedNotifications<std::vector<PublicationFieldChange>>> notifications;
};

/* A detached candidate may prepare geometry with the node's ordinary builder.
 * Publish its scalar state and complete child list without duplicating the
 * commit, callback-drain and exception-ordering protocol at every node API. */
inline void
publish_scalar_fields_and_children(SoSeparator &target,
	const SoSeparator &candidate)
{
    std::vector<SoNode *> children;
    children.reserve(static_cast<size_t>(candidate.getNumChildren()));
    for (int i = 0; i < candidate.getNumChildren(); ++i)
	children.push_back(candidate.getChild(i));

    auto replacement = target.getChildren()->prepareReplacement(children);
    PreparedScalarFields fields(target, candidate);
    replacement->commit();
    fields.commit();

    fields.restore();
    std::exception_ptr failure;
    try { replacement->notify(); }
    catch (...) { failure = std::current_exception(); }
    if (fields.changed())
	fields.notify(failure);
    if (failure)
	std::rethrow_exception(failure);
}

class PreparedScalarNodes {
public:
    void reserve(size_t count) { this->nodes.reserve(count); }
    void prepare(SoNode &target, const SoNode &candidate)
    {
	this->nodes.push_back(
	    std::make_unique<PreparedScalarFields>(target, candidate));
    }
    void commit()
    {
	for (auto &node : this->nodes)
	    node->commit();
    }
    void restore()
    {
	for (auto &node : this->nodes)
	    node->restore();
    }
    void notify(std::exception_ptr &failure)
    {
	for (auto &node : this->nodes)
	    node->notify(failure);
    }
private:
    std::vector<std::unique_ptr<PreparedScalarFields>> nodes;
};

#endif /* BOBOL_SCALAR_PUBLICATION_PRIVATE_H */
