/*              C A D _ P U B L I C A T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "cad_publication_private.h"
#include "transaction_fault_private.h"

#include "bu/log.h"

#include <Obol/cad/CadSceneValidation.h>
#include <Obol/cad/CadSceneReplacement.h>
#include <Obol/cad/SoCADAssembly.h>

#include <utility>

namespace {

const char *
publication_operation(const char *operation)
{
    return operation ? operation : "unknown";
}

bool
report_geometry_validation(const Obol::CadGeometryValidation &validation,
    const char *operation)
{
    if (validation)
	return true;
    bu_log("libBObol: rejected %s CAD geometry: %s at update %zu, "
	"element %zu\n", publication_operation(operation),
	Obol::cadGeometryErrorName(validation.error), validation.updateIndex,
	validation.elementIndex);
    return false;
}

bool
report_scene_validation(const Obol::CadSceneValidation &validation,
    const char *operation)
{
    if (validation)
	return true;
    bu_log("libBObol: rejected %s CAD scene mutation: %s at update %zu\n",
	publication_operation(operation),
	Obol::cadSceneErrorName(validation.error), validation.updateIndex);
    return false;
}

bool
report_scene_replacement(const Obol::CadSceneReplacementResult &result,
    const char *operation)
{
    if (result)
	return true;
    switch (result.error) {
	case Obol::CadSceneReplacementError::Geometry:
	    return report_geometry_validation(result.geometry, operation);
	case Obol::CadSceneReplacementError::Instances:
	    return report_scene_validation(result.instances, operation);
	case Obol::CadSceneReplacementError::ResourceUnavailable:
	    break;
	case Obol::CadSceneReplacementError::Valid:
	    break;
    }
    bu_log("libBObol: rejected %s CAD scene replacement: %s\n",
	publication_operation(operation),
	Obol::cadSceneReplacementErrorName(result.error));
    return false;
}

bool
report_scene_mutation(const Obol::CadSceneMutationResult &result,
	const Obol::CadSceneMutation &mutation, const char *operation)
{
    if (result)
	return true;
    if (result.domain == Obol::CadSceneMutationDomain::Parts)
	return report_geometry_validation(result.geometry, operation);
    bu_log("libBObol: rejected %s CAD scene transaction: %s/%s "
	"at update %zu\n", publication_operation(operation),
	Obol::cadSceneMutationDomainName(result.domain),
	Obol::cadSceneErrorName(result.scene.error), result.scene.updateIndex);
    if (result.domain == Obol::CadSceneMutationDomain::Cuts &&
	result.scene.updateIndex < mutation.cuts.size()) {
	const Obol::InstanceId instance =
	    mutation.cuts[result.scene.updateIndex].instance;
	bu_log("libBObol: rejected cut target instance=%016llx%016llx\n",
	    static_cast<unsigned long long>(instance.w0),
	    static_cast<unsigned long long>(instance.w1));
    } else if (result.domain == Obol::CadSceneMutationDomain::Styles &&
	result.scene.updateIndex < mutation.styles.size()) {
	const Obol::InstanceId instance =
	    mutation.styles[result.scene.updateIndex].instance;
	bu_log("libBObol: rejected style target instance=%016llx%016llx\n",
	    static_cast<unsigned long long>(instance.w0),
	    static_cast<unsigned long long>(instance.w1));
    }
    return false;
}

std::shared_ptr<const Obol::PartGeometry>
build_geometry(Obol::PartGeometryBuilder geometry, const char *operation,
    Obol::CadAggregateProxyPolicy proxyPolicy)
{
    Obol::CadGeometryAdmission admission = Obol::cadAdmitPartGeometry(
	std::move(geometry), proxyPolicy);
    if (!report_geometry_validation(admission.validation, operation))
	return std::shared_ptr<const Obol::PartGeometry>();
    return admission.geometry.shared();
}

} // namespace

bool
bobol_cad_admit_geometry(
    const std::shared_ptr<const Obol::PartGeometry> &geometry,
    Obol::ValidatedPartGeometry &admitted, const char *operation)
{
    if (!geometry) {
	Obol::CadGeometryValidation validation;
	validation.error = Obol::CadGeometryError::NullGeometry;
	report_geometry_validation(validation, operation);
	return false;
    }
    admitted = Obol::ValidatedPartGeometry(geometry);
    return true;
}

std::shared_ptr<const Obol::PartGeometry>
bobol_cad_build_geometry(Obol::PartGeometryBuilder geometry,
    const char *operation)
{
    return build_geometry(std::move(geometry), operation,
	Obol::CadAggregateProxyPolicy::Strict);
}

std::shared_ptr<const Obol::PartGeometry>
bobol_cad_build_geometry_with_optional_proxy(
    Obol::PartGeometryBuilder geometry, const char *operation)
{
    return build_geometry(std::move(geometry), operation,
	Obol::CadAggregateProxyPolicy::DiscardInvalid);
}

bool
bobol_cad_validate_shared_parts(
    const std::vector<Obol::PartUpdate> &updates,
    const char *operation)
{
    return report_geometry_validation(
	Obol::cadValidatePartUpdates(updates), operation);
}

bool
bobol_cad_validate_instances(
    const std::vector<Obol::InstanceUpdate> &updates,
    const char *operation)
{
    return report_scene_validation(
	Obol::cadValidateInstanceUpdates(updates), operation);
}

bool
bobol_cad_validate_instance(const Obol::InstanceUpdate &update,
    size_t updateIndex, const char *operation)
{
    Obol::CadSceneValidation validation =
	Obol::cadValidateInstanceRecord(update.instance, update.record);
    if (validation)
	return true;
    validation.updateIndex = updateIndex;
    return report_scene_validation(validation, operation);
}

bool
bobol_cad_validate_styles(
    const std::vector<Obol::InstanceStyleUpdate> &updates,
    const char *operation)
{
    return report_scene_validation(
	Obol::cadValidateInstanceStyleUpdates(updates), operation);
}

bool
bobol_cad_publish_styles(SoCADAssembly *assembly,
    const std::vector<Obol::InstanceStyleUpdate> &updates,
    const char *operation)
{
    if (!assembly) {
	bu_log("libBObol: rejected %s CAD scene mutation: null assembly\n",
	    publication_operation(operation));
	return false;
    }
    return report_scene_validation(
	assembly->updateInstanceStyles(updates), operation);
}

bool
bobol_cad_publish_mutation(SoCADAssembly *assembly,
    const Obol::CadSceneMutation &mutation, const char *operation)
{
    if (!assembly) {
	bu_log("libBObol: rejected %s CAD scene transaction: null assembly\n",
	    publication_operation(operation));
	return false;
    }
    if (bobol_transaction_fault_requested(
	    BObolTransactionFaultPoint::RETAINED_SCENE_COMMIT)) {
	bu_log("libBObol: injected resource denial before %s CAD scene "
	    "transaction commit\n", publication_operation(operation));
	return false;
    }
    return report_scene_mutation(
	assembly->applySceneMutation(mutation), mutation, operation);
}

bool
bobol_cad_validate_mutation(const SoCADAssembly *assembly,
    const Obol::CadSceneMutation &mutation, const char *operation)
{
    if (!assembly) {
	bu_log("libBObol: rejected %s CAD scene transaction: null assembly\n",
	    publication_operation(operation));
	return false;
    }
    return report_scene_mutation(
	assembly->validateSceneMutation(mutation), mutation, operation);
}

bool
bobol_cad_replace_scene(SoCADAssembly *assembly,
    const std::vector<Obol::PartUpdate> &parts,
    const std::vector<Obol::InstanceUpdate> &instances,
    const char *operation)
{
    if (!assembly) {
	bu_log("libBObol: rejected %s CAD scene replacement: null assembly\n",
	    publication_operation(operation));
	return false;
    }
    if (bobol_transaction_fault_requested(
	    BObolTransactionFaultPoint::RETAINED_SCENE_COMMIT)) {
	bu_log("libBObol: injected resource denial before %s CAD scene "
	    "replacement commit\n", publication_operation(operation));
	return false;
    }
    return report_scene_replacement(
	assembly->replaceScene(parts, instances), operation);
}
