/*              T E S T _ C O M P A C T _ E D I T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol.h"
#include "BObol/BEditPreview.h"
#include "BObol/BMeshShape.h"
#include "BObol/BVListShape.h"
#include "Obol/cad/SoCADAssembly.h"
#include "bu/app.h"
#include "bu/str.h"

#include <Inventor/nodes/SoGroup.h>

#include <memory>
#include <stdio.h>
#include <vector>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

static SbBool
near_point(const SbVec3f &a, const SbVec3f &b)
{
    return (a - b).length() < 0.0001f ? TRUE : FALSE;
}

static SbBool
same_compact_handle(const BObolCompactInstanceHandle &a,
	const BObolCompactInstanceHandle &b)
{
    return a.sourceNodeId == b.sourceNodeId &&
	a.instanceWord0 == b.instanceWord0 &&
	a.instanceWord1 == b.instanceWord1 ? TRUE : FALSE;
}

static int
count_nodes_of_type(SoNode *node, SoType type)
{
    if (!node)
	return 0;
    int count = node->isOfType(type) ? 1 : 0;
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return count;
    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	count += count_nodes_of_type(group->getChild(i), type);
    return count;
}

static SoBRLVListShape *
preview_line_shape(SoNode *node)
{
    if (!node)
	return NULL;
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<SoBRLVListShape *>(node);
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return NULL;
    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoBRLVListShape *shape = preview_line_shape(group->getChild(i));
	if (shape)
	    return shape;
    }
    return NULL;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    if (argc != 1)
	FAIL("unexpected arguments");
    bobol_init(NULL);

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->ref();
    source->path = "assembly.g";
    source->selectedColor = SbColor(0.9f, 0.2f, 0.1f);
    if (!source->setDisplayState(TRUE, 17, 23, TRUE, FALSE, FALSE,
	    1, 4, 0.25f, TRUE, SbColor(0.1f, 0.3f, 0.7f),
	    FALSE, SbColor(1.0f, 1.0f, 1.0f), 5)) {
	source->unref();
	FAIL("source display setup should change state");
    }

    SbMatrix sourcePlacement;
    sourcePlacement.setTranslate(SbVec3f(10.0f, 20.0f, 30.0f));
    if (!source->setPlacementState(TRUE, sourcePlacement, FALSE,
	    SbVec3f(0.0f, 0.0f, 0.0f), FALSE, 0.0f)) {
	source->unref();
	FAIL("source placement setup should change state");
    }

    std::shared_ptr<Obol::PartGeometry> geometry =
	std::make_shared<Obol::PartGeometry>();
    Obol::WireRep wire;
    wire.segmentPoints.push_back(SbVec3f(1.0f, 2.0f, 3.0f));
    wire.segmentPoints.push_back(SbVec3f(4.0f, 2.0f, 3.0f));
    wire.segmentIds.push_back(11);
    wire.bounds.makeEmpty();
    wire.bounds.extendBy(wire.segmentPoints[0]);
    wire.bounds.extendBy(wire.segmentPoints[1]);
    geometry->wire = wire;

    BObolCompactOccurrence occurrence;
    occurrence.geometry = geometry;
    occurrence.summary.valid = TRUE;
    occurrence.summary.shapeKind = BObolRealizedShapeSummary::SHAPE_VLIST;
    occurrence.summary.path = "assembly.g/part.s";
    occurrence.summary.sourceName = "part.s";
    occurrence.summary.sourceType = "wire";
    occurrence.summary.visible = TRUE;
    occurrence.summary.selectable = TRUE;
    occurrence.localTransform.setTranslate(SbVec3f(2.0f, 3.0f, 4.0f));
    occurrence.lodBacked = TRUE;
    occurrence.occurrenceIndex = 3;

    /* Asynchronous draw roots can be attached after the semantic selection
     * delta was emitted.  They retain the current path rule while empty, and
     * every subsequently streamed occurrence must inherit it immediately. */
    std::vector<SbString> selectedPaths = {
	SbString("assembly.g/part.s")
    };
    if (source->syncCompactInstanceSelectedPaths(selectedPaths) <= 0) {
	source->unref();
	FAIL("empty compact source should retain semantic selection paths");
    }
    if (source->setCompactOccurrence(occurrence) != 1) {
	source->unref();
	FAIL("compact occurrence setup should succeed");
    }

    BObolCompactInstanceHandle pathHandle;
    BObolCompactInstanceSummary pathSummary;
    if (!source->getCompactInstanceForPath("assembly.g/part.s", FALSE,
	    TRUE, pathHandle, pathSummary) || !pathSummary.valid ||
	!BU_STR_EQUAL(pathSummary.path.getString(), "assembly.g/part.s") ||
	source->getCompactInstanceForPath("assembly.g/missing.s", FALSE,
	    TRUE, pathHandle, pathSummary)) {
	source->unref();
	FAIL("compact semantic paths should resolve without materialization");
    }

    SoGroup *sceneRoot = new SoGroup;
    sceneRoot->ref();
    sceneRoot->addChild(source);
    BObolSceneController scene(sceneRoot);
    BObolSceneTreeSummary treeSummary;
    if (!scene.getSceneTreeSummaryForPath("assembly.g/part.s",
	    treeSummary) || !treeSummary.valid ||
	!treeSummary.isDatabaseSource ||
	!BU_STR_EQUAL(treeSummary.path.getString(), "assembly.g/part.s") ||
	!BU_STR_EQUAL(treeSummary.ownerSourcePath.getString(), "assembly.g") ||
	treeSummary.ownerSourceInstanceKey.getLength() == 0) {
	sceneRoot->unref();
	source->unref();
	FAIL("scene paths should address compact occurrences without leaf nodes");
    }

    if (count_nodes_of_type(source, SoBRLVListShape::getClassTypeId()) != 0 ||
	count_nodes_of_type(source, SoBRLMeshShape::getClassTypeId()) != 0) {
	source->unref();
	FAIL("normal compact state should retain no off-scene render shape");
    }

    BObolCompactInstanceHandle compact;
    if (!source->getCompactInstanceHandle(0, compact)) {
	source->unref();
	FAIL("streamed selected compact occurrence should have a handle");
    }

    BObolCompactInstanceSummary visibilitySummary;
    if (source->setCompactInstanceVisibilityOverrideForPathMatch(
	    "/assembly.g/part.s", BOBOL_COMPACT_PATH_EXACT, FALSE) <= 0 ||
	!source->getCompactInstanceSummary(compact, visibilitySummary) ||
	visibilitySummary.visible ||
	source->setCompactInstanceVisibilityOverrideForPathMatch(
	    "assembly.g/part.s", BOBOL_COMPACT_PATH_EXACT, TRUE) <= 0 ||
	!source->getCompactInstanceSummary(compact, visibilitySummary) ||
	!visibilitySummary.visible) {
	sceneRoot->unref();
	source->unref();
	FAIL("compact presentation paths should normalize leading slashes");
    }

    BObolCompactInstanceSummary before;
    if (!source->getCompactInstanceSummary(compact, before) ||
	!before.selected || !before.lodBacked ||
	!before.appearanceColorValid ||
	!near_point(before.appearanceColor,
	    SbColor(0.9f, 0.2f, 0.1f)) ||
	before.lineStyle != 1 || before.lineWidth != 4 ||
	!NEAR_EQUAL(before.transparency, 0.25f, SMALL_FASTF)) {
	source->unref();
	FAIL("compact summary should expose effective selected appearance and LoD");
    }

    SbVec3f expectedStart;
    before.localToSource.multVecMatrix(wire.segmentPoints[0], expectedStart);
    if (near_point(expectedStart, wire.segmentPoints[0])) {
	source->unref();
	FAIL("compact local-to-source summary should include source placement");
    }

    BObolFeatureOwner owner;
    owner.ownerToken = source;
    owner.ownerId = "compact-edit-test";
    owner.ownerRole = "edit";
    BObolFeatureStore features;
    BObolFeatureHandle preview = features.promoteCompactInstanceForEdit(
	"part-edit", *source, compact, "move-part", "transform", &owner);
    if (!preview.isValid()) {
	source->unref();
	FAIL("compact promotion should publish an edit preview");
    }

    BObolCompactInstanceHandle boundHandle;
    BObolCompactInstanceSummary bound;
    if (!features.compactEditBinding(preview, boundHandle, bound) ||
	!same_compact_handle(boundHandle, compact) ||
	bound.geometryIdentity != before.geometryIdentity ||
	bound.geometryRevision != before.geometryRevision ||
	bound.appearanceRevision != before.appearanceRevision ||
	bound.placementRevision != before.placementRevision ||
	bound.selectionRevision != before.selectionRevision ||
	!bound.localToSource.equals(before.localToSource, 0.000001f) ||
	!bound.lodBacked || !bound.selected) {
	source->unref();
	FAIL("edit binding should retain compact identity, revisions, placement, selection, and LoD");
    }

    std::vector<SbVec3f> previewPoints;
    BObolFeatureStyle previewStyle;
    if (!features.points(preview, previewPoints) ||
	previewPoints.size() != 2 ||
	!near_point(previewPoints[0], expectedStart) ||
	!features.style(preview, previewStyle) ||
	!previewStyle.hasColor ||
	!near_point(previewStyle.color, before.appearanceColor) ||
	!previewStyle.hasLineWidth || previewStyle.lineWidth != 4 ||
	!previewStyle.hasLineStyle || previewStyle.lineStyle != 1 ||
	!previewStyle.hasTransparency ||
	!NEAR_EQUAL(previewStyle.transparency, 0.25f, SMALL_FASTF)) {
	source->unref();
	FAIL("promoted preview should preserve transformed geometry and appearance");
    }

    SoNode *previewNode = features.node(preview);
    SoBRLVListShape *previewShape = preview_line_shape(previewNode);
    if (!previewNode ||
	!previewNode->isOfType(SoBRLEditPreview::getClassTypeId()) ||
	!previewShape || previewShape->lineWidth.getValue() != 4 ||
	previewShape->lineStyle.getValue() != 1 ||
	!NEAR_EQUAL(previewShape->transparency.getValue(), 0.25f,
	    SMALL_FASTF)) {
	source->unref();
	FAIL("promoted Coin preview should receive the effective compact style");
    }

    BObolCompactInstanceSummary during;
    if (!source->getCompactInstanceSummary(compact, during) ||
	during.geometryIdentity != before.geometryIdentity ||
	during.placementRevision != before.placementRevision ||
	count_nodes_of_type(source, SoBRLVListShape::getClassTypeId()) != 0) {
	source->unref();
	FAIL("promotion should not mutate or materialize the compact source");
    }

    BObolCompactInstanceHandle wrong = compact;
    wrong.instanceWord1 ^= 1u;
    if (features.demoteCompactInstanceFromEdit(preview, wrong) ||
	!features.existsOwned("part-edit", BOBOL_FEATURE_SCOPE_LOCAL, &owner)) {
	source->unref();
	FAIL("demotion should reject a different compact identity");
    }
    if (!features.demoteCompactInstanceFromEdit(preview, compact) ||
	features.existsOwned("part-edit", BOBOL_FEATURE_SCOPE_LOCAL, &owner)) {
	source->unref();
	FAIL("demotion should remove only the transient edit preview");
    }

    BObolCompactInstanceSummary after;
    if (!source->getCompactInstanceSummary(compact, after) ||
	after.geometryIdentity != before.geometryIdentity ||
	after.appearanceRevision != before.appearanceRevision ||
	after.placementRevision != before.placementRevision ||
	after.selectionRevision != before.selectionRevision ||
	!after.localToSource.equals(before.localToSource, 0.000001f) ||
	!after.lodBacked || !after.selected ||
	count_nodes_of_type(source, SoBRLVListShape::getClassTypeId()) != 0) {
	source->unref();
	FAIL("demotion should retain compact identity, appearance, placement, selection, and LoD");
    }

    sceneRoot->unref();
    source->unref();
    return 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
