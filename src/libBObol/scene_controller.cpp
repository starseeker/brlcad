/*                S C E N E _ C O N T R O L L E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include <Inventor/tools/SbModernUtils.h>

#include "bu/str.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BDrawCache.h"
#include "BObol/BGrid.h"
#include "BObol/BLineLayerOverlay.h"
#include "BObol/BMaterialObject.h"
#include "BObol/BMeshShape.h"
#include "BObol/BRealizeAction.h"
#include "BObol/BSceneController.h"
#include "BObol/BSceneGroup.h"
#include "BObol/BVListShape.h"
#include "identity_counter_private.h"
#include "performance_private.h"
#include "scalar_publication_private.h"
#include "database_source_private.h"
#include "database_source_realization.h"

#include "raytrace.h"

#include <Inventor/SbName.h>
#include <Inventor/misc/SoChildList.h>
#include <exception>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/lists/SoAuditorList.h>
#include <Inventor/lists/SbList.h>
#include <Inventor/misc/SoNotRec.h>
#include <Inventor/sensors/SoNodeSensor.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoMatrixTransform.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSeparator.h>

#include <algorithm>
#include <limits>
#include <set>
#include <stdexcept>
#include <map>
#include <string.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <type_traits>
#include <utility>
#include <vector>

static const char *
skip_leading_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

static std::string
normalized_index_key(const char *key)
{
    const char *normalized = skip_leading_slash(key);
    const size_t len = strlen(normalized);
    std::vector<char> stable;
    stable.resize(len + 1);
    if (len)
	memcpy(stable.data(), normalized, len);
    stable[len] = '\0';
    return std::string(stable.data(), len);
}

static int
scene_path_equal(const char *stored, const char *path)
{
    if (!stored || !path)
	return 0;
    if (bu_strcmp(stored, path) == 0)
	return 1;
    return bu_strcmp(skip_leading_slash(stored), skip_leading_slash(path)) == 0;
}

static bool
scene_path_contains(const std::string &ancestor, const std::string &path)
{
    return ancestor == path ||
	(path.size() > ancestor.size() &&
	 path.compare(0, ancestor.size(), ancestor) == 0 &&
	 path[ancestor.size()] == '/');
}

static int
database_source_path_equal(const SoBRLDatabaseSource *source,
			   const char *path)
{
    if (!source || !path)
	return 0;
    const char *sourcePath = source->path.getValue().getString();
    return scene_path_equal(sourcePath, path);
}

static SbString
database_source_effective_instance_key(const SoBRLDatabaseSource *source)
{
    if (!source)
	return "";
    const SbString key = source->instanceKey.getValue();
    if (key.getLength() > 0)
	return key;
    return source->path.getValue();
}

static int
database_source_instance_key_equal(const SoBRLDatabaseSource *source,
				   const char *instanceKey)
{
    if (!source || !instanceKey)
	return 0;
    const SbString key = database_source_effective_instance_key(source);
    const char *stored = key.getString();
    return scene_path_equal(stored, instanceKey);
}

static void
index_put(std::unordered_map<std::string, SoBRLDatabaseSource *> &index,
	  const char *key,
	  SoBRLDatabaseSource *source)
{
    if (!key || !key[0] || !source)
	return;

    const char *normalized = skip_leading_slash(key);
    if (normalized && normalized[0])
	index[std::string(normalized)] = source;
}

static void
parent_index_put(std::unordered_map<std::string, SoGroup *> &index,
		 const char *key,
		 SoGroup *parent)
{
    if (!key || !key[0] || !parent)
	return;

    const char *normalized = skip_leading_slash(key);
    if (normalized && normalized[0])
	index[std::string(normalized)] = parent;
}

static void
group_index_put(std::unordered_map<std::string, SoGroup *> &index,
		const char *key,
		SoGroup *group)
{
    if (!key || !key[0] || !group)
	return;

    const char *normalized = skip_leading_slash(key);
    if (normalized && normalized[0])
	index[std::string(normalized)] = group;
}

static SbString
scene_group_index_path(const SoGroup *group)
{
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return "";

    const SoBRLSceneGroup *sceneGroup =
	static_cast<const SoBRLSceneGroup *>(group);
    return sceneGroup->groupPath.getValue();
}

static void
source_order_put(std::vector<SoBRLDatabaseSource *> &order,
		 std::unordered_map<SoBRLDatabaseSource *, size_t> &orderIndex,
		 SoBRLDatabaseSource *source)
{
    if (!source || orderIndex.find(source) != orderIndex.end())
	return;

    orderIndex[source] = order.size();
    order.push_back(source);
}

static SoBRLDatabaseSource *
find_database_source_instance_recursive(SoGroup *group,
					const char *sourceInstanceKey,
					SoGroup **parentOut = NULL, int *childIndexOut = NULL)
{
    if (parentOut)
	*parentOut = NULL;
    if (childIndexOut)
	*childIndexOut = -1;
    if (!group || !sourceInstanceKey)
	return NULL;

    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node)
	    continue;
	if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    SoBRLDatabaseSource *source =
		static_cast<SoBRLDatabaseSource *>(node);
	    if (database_source_instance_key_equal(source,
						   sourceInstanceKey)) {
		if (parentOut)
		    *parentOut = group;
		if (childIndexOut)
		    *childIndexOut = i;
		return source;
	    }
	}
	if (node->isOfType(SoGroup::getClassTypeId())) {
	    SoBRLDatabaseSource *found =
		find_database_source_instance_recursive(
		    static_cast<SoGroup *>(node), sourceInstanceKey,
		    parentOut, childIndexOut);
	    if (found)
		return found;
	}
    }

    return NULL;
}

/* Small parent chains use the native list's inline storage. Once a traversal
 * exceeds that capacity, hashed membership preserves bounded graph-walk cost. */
class SceneNodeSet : private SbList<SoNode *> {
public:
    bool insert(SoNode *node)
    {
	if (this->expanded) return this->nodes.insert(node).second;
	if (this->find(node) >= 0) return false;
	if (this->getLength() < this->getArraySize()) {
	    this->append(node);
	    return true;
	}
	this->nodes.reserve(size_t(this->getLength()) + 1);
	for (int i = 0; i < this->getLength(); ++i) this->nodes.insert((*this)[i]);
	this->nodes.insert(node);
	this->truncate(0);
	this->expanded = true;
	return true;
    }
    void erase(SoNode *node)
    {
	if (this->expanded) this->nodes.erase(node);
	else {
	    const int index = this->find(node);
	    if (index >= 0) this->removeFast(index);
	}
    }
private:
    std::unordered_set<SoNode *> nodes;
    bool expanded = false;
};

template <typename Visit, typename Children>
static bool
scene_walk(SoNode *root, SoGroup *parent, Visit visit, Children children)
{
    struct Entry { SoNode *node; SoGroup *parent; bool exit; };
    SbList<Entry> pending;
    pending.append({root, parent, false});
    SceneNodeSet ancestors;
    while (pending.getLength() > 0) {
        const auto entry = pending.pop();
        if (!entry.node) continue;
        if (entry.exit) { ancestors.erase(static_cast<SoGroup *>(entry.node)); continue; }
        visit(entry.node, entry.parent);
        if (!entry.node->isOfType(SoGroup::getClassTypeId())) continue;
        auto *group = static_cast<SoGroup *>(entry.node);
        if (!ancestors.insert(group)) return false;
        pending.append({group, entry.parent, true});
        const auto *prepared = children(group);
        if (prepared) {
	for (auto it = prepared->rbegin(); it != prepared->rend(); ++it) pending.append({*it, group, false});
        } else {
	for (int i = group->getNumChildren(); i > 0; --i) pending.append({group->getChild(i - 1), group, false});
        }
    }
    return true;
}

static void
index_database_sources_recursive(
    SoGroup *group,
    std::unordered_map<std::string, SoGroup *> &groupPathIndex,
    std::unordered_map<std::string, SoBRLDatabaseSource *> &pathIndex,
    std::unordered_map<std::string, std::vector<SoBRLDatabaseSource *> > &pathInstancesIndex,
    std::unordered_map<std::string, SoBRLDatabaseSource *> &instanceIndex,
    std::unordered_map<std::string, SoGroup *> &instanceParentIndex,
    std::vector<SoBRLDatabaseSource *> &sourceOrder,
    std::unordered_map<SoBRLDatabaseSource *, size_t> &sourceOrderIndex)
{
    if (!scene_walk(group, nullptr, [&](SoNode *node, SoGroup *parent) {
	if (node->isOfType(SoGroup::getClassTypeId())) {
	    auto *nodeGroup = static_cast<SoGroup *>(node);
	    const auto path = scene_group_index_path(nodeGroup);
	    group_index_put(groupPathIndex, path.getString(), nodeGroup);
	}
	if (!parent || !node->isOfType(SoBRLDatabaseSource::getClassTypeId())) return;
	auto *source = static_cast<SoBRLDatabaseSource *>(node);
	const auto path = source->path.getValue();
	const auto key = database_source_effective_instance_key(source);
	index_put(pathIndex, path.getString(), source);
	if (!sourceOrderIndex.count(source)) pathInstancesIndex[normalized_index_key(path.getString())].push_back(source);
	index_put(instanceIndex, key.getString(), source);
	parent_index_put(instanceParentIndex, key.getString(), parent);
	source_order_put(sourceOrder, sourceOrderIndex, source);
    }, [](SoGroup *) -> const std::vector<SoNode *> * { return nullptr; }))
	throw std::invalid_argument("scene root contains a cycle");
}

static int
count_database_sources_recursive(const SoGroup *group)
{
    if (!group)
	return 0;

    int count = 0;
    for (int i = 0; i < group->getNumChildren(); i++) {
	const SoNode *node = group->getChild(i);
	if (!node)
	    continue;
	if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    count++;
	}
	if (node->isOfType(SoGroup::getClassTypeId()))
	    count += count_database_sources_recursive(
			 static_cast<const SoGroup *>(node));
    }
    return count;
}

static int
count_scene_groups_recursive(const SoGroup *group)
{
    if (!group)
	return 0;

    int count = 0;
    for (int i = 0; i < group->getNumChildren(); i++) {
	const SoNode *node = group->getChild(i);
	if (!node)
	    continue;
	if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	    continue;
	if (node->isOfType(SoBRLSceneGroup::getClassTypeId()))
	    count++;
	if (node->isOfType(SoGroup::getClassTypeId()))
	    count += count_scene_groups_recursive(
			 static_cast<const SoGroup *>(node));
    }
    return count;
}

static void
scene_group_path_components(const char *groupPath,
			    std::vector<std::string> &components)
{
    components.clear();
    if (!groupPath)
	return;

    const char *cursor = groupPath;
    while (*cursor == '/')
	cursor++;

    while (*cursor) {
	const char *start = cursor;
	while (*cursor && *cursor != '/')
	    cursor++;
	if (cursor > start)
	    components.push_back(std::string(start,
					     static_cast<size_t>(cursor - start)));
	while (*cursor == '/')
	    cursor++;
    }
}

static SoGroup *
scene_root_group(SoNode *node)
{
    if (!node || !node->isOfType(SoGroup::getClassTypeId()))
	return NULL;
    return static_cast<SoGroup *>(node);
}

static const SoGroup *
scene_root_group_const(const SoNode *node)
{
    if (!node || !node->isOfType(SoGroup::getClassTypeId()))
	return NULL;
    return static_cast<const SoGroup *>(node);
}

static int
scene_group_node_name_equal(const SoNode *node, const char *leafName)
{
    if (!node || !leafName)
	return 0;
    const SbName nodeName = node->getName();
    if (bu_strcmp(nodeName.getString(), leafName) == 0)
	return 1;

    if (node->isOfType(SoBRLSceneGroup::getClassTypeId())) {
	const SoBRLSceneGroup *group =
	    static_cast<const SoBRLSceneGroup *>(node);
	const char *groupPath = group->groupPath.getValue().getString();
	const char *groupLeaf = strrchr(groupPath, '/');
	groupLeaf = groupLeaf ? groupLeaf + 1 : groupPath;
	if (groupLeaf && bu_strcmp(groupLeaf, leafName) == 0)
	    return 1;
    }

    return 0;
}

static SoGroup *
scene_group_find_child(SoGroup *parent, const char *leafName)
{
    if (!parent || !leafName || !leafName[0])
	return NULL;

    for (int i = 0; i < parent->getNumChildren(); i++) {
	SoNode *child = parent->getChild(i);
	if (!child || !child->isOfType(SoBRLSceneGroup::getClassTypeId()) ||
	    !scene_group_node_name_equal(child, leafName))
	    continue;
	return static_cast<SoGroup *>(child);
    }

    return NULL;
}

static int
scene_group_find_child_index(SoGroup *parent, const SoNode *child)
{
    if (!parent || !child)
	return -1;
    for (int i = 0; i < parent->getNumChildren(); i++) {
	if (parent->getChild(i) == child)
	    return i;
    }
    return -1;
}

static const SoGroup *
scene_group_find_child_const(const SoGroup *parent, const char *leafName)
{
    if (!parent || !leafName || !leafName[0])
	return NULL;

    for (int i = 0; i < parent->getNumChildren(); i++) {
	const SoNode *child = parent->getChild(i);
	if (!child || !child->isOfType(SoBRLSceneGroup::getClassTypeId()) ||
	    !scene_group_node_name_equal(child, leafName))
	    continue;
	return static_cast<const SoGroup *>(child);
    }

    return NULL;
}

static const SoGroup *
scene_group_find_path_const(const SoNode *sceneRoot, const char *groupPath)
{
    if (!groupPath)
	return NULL;

    const SoGroup *current = scene_root_group_const(sceneRoot);
    if (!current)
	return NULL;

    std::vector<std::string> components;
    scene_group_path_components(groupPath, components);
    for (size_t i = 0; i < components.size(); i++) {
	current = scene_group_find_child_const(current,
					       components[i].c_str());
	if (!current)
	    return NULL;
    }

    return current;
}

static SbString
scene_child_summary_path(const SbString &parentPath, const SoNode *child)
{
    if (!child)
	return "";

    const SbName childName = child->getName();
    if (childName.getLength() == 0)
	return "";

    if (parentPath.getLength() == 0 ||
	bu_strcmp(parentPath.getString(), "/") == 0)
	return SbString(childName.getString());

    SbString childPath = parentPath;
    childPath += "/";
    childPath += childName.getString();
    return childPath;
}

static SbString
scene_group_append_path(const SbString &parentPath, const char *leafName)
{
    if (!leafName || !leafName[0])
	return parentPath;

    if (parentPath.getLength() == 0 ||
	bu_strcmp(parentPath.getString(), "/") == 0)
	return SbString(leafName);

    SbString childPath = parentPath;
    childPath += "/";
    childPath += leafName;
    return childPath;
}

static SbString
scene_group_component_path(const std::vector<std::string> &components,
			   size_t count)
{
    SbString path("");
    const size_t limit = count < components.size() ?
			 count : components.size();
    for (size_t i = 0; i < limit; i++)
	path = scene_group_append_path(path, components[i].c_str());
    return path;
}

static SbString
scene_group_summary_path(const SoNode *node, const SbString &fallbackPath)
{
    if (node && node->isOfType(SoBRLSceneGroup::getClassTypeId())) {
	const SoBRLSceneGroup *group =
	    static_cast<const SoBRLSceneGroup *>(node);
	if (group->groupPath.getValue().getLength() > 0)
	    return group->groupPath.getValue();
    }

    if (node && node->isOfType(SoBRLLineLayerOverlay::getClassTypeId())) {
	const SoBRLLineLayerOverlay *overlay =
	    static_cast<const SoBRLLineLayerOverlay *>(node);
	if (overlay->overlayId.getValue().getLength() > 0)
	    return overlay->overlayId.getValue();
    }

    if (node && node->isOfType(SoBRLGrid::getClassTypeId())) {
	const SoBRLGrid *grid = static_cast<const SoBRLGrid *>(node);
	if (grid->overlayId.getValue().getLength() > 0)
	    return grid->overlayId.getValue();
    }

    return fallbackPath;
}

static constexpr float scene_state_tolerance = 0.000001f;

static int
scene_group_float_different(float a, float b)
{
    const float delta = a - b;
    return delta < -scene_state_tolerance || delta > scene_state_tolerance;
}

static int
scene_group_color_equal(const SbColor &a, const SbColor &b)
{
    return !scene_group_float_different(a[0], b[0]) &&
	   !scene_group_float_different(a[1], b[1]) &&
	   !scene_group_float_different(a[2], b[2]);
}

static int
scene_group_vec3f_equal(const SbVec3f &a, const SbVec3f &b)
{
    return !scene_group_float_different(a[0], b[0]) &&
	   !scene_group_float_different(a[1], b[1]) &&
	   !scene_group_float_different(a[2], b[2]);
}

static int
scene_node_is_shape(const SoNode *node)
{
    return node &&
	   (node->isOfType(SoBRLVListShape::getClassTypeId()) ||
	    node->isOfType(SoBRLMeshShape::getClassTypeId()));
}

static SbString
scene_shape_node_path(const SoNode *node)
{
    if (!node)
	return "";
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<const SoBRLVListShape *>(node)->sourcePath.getValue();
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<const SoBRLMeshShape *>(node)->sourcePath.getValue();
    return "";
}

static int
scene_shape_path_equal(const SoNode *node, const char *shapePath)
{
    if (!node || !shapePath || !shapePath[0])
	return 0;

    const SbString nodePath = scene_shape_node_path(node);
    if (nodePath.getLength() == 0)
	return 0;

    const char *stored = nodePath.getString();
    if (bu_strcmp(stored, shapePath) == 0)
	return 1;
    return bu_strcmp(skip_leading_slash(stored),
		  skip_leading_slash(shapePath)) == 0;
}

static SoNode *
scene_shape_find_in_group(SoGroup *group, const char *shapePath,
			  SoGroup **parentOut)
{
    if (!group || !shapePath)
	return NULL;

    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *child = group->getChild(i);
	if (!child)
	    continue;

	if (scene_node_is_shape(child) &&
	    scene_shape_path_equal(child, shapePath)) {
	    if (parentOut)
		*parentOut = group;
	    return child;
	}

	if (!child->isOfType(SoGroup::getClassTypeId()) ||
	    child->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	    continue;

	SoNode *found = scene_shape_find_in_group(
			    static_cast<SoGroup *>(child), shapePath, parentOut);
	if (found)
	    return found;
    }

    return NULL;
}

static SoNode *
scene_shape_find_path(SoNode *sceneRoot, const char *shapePath,
		      SoGroup **parentOut)
{
    if (parentOut)
	*parentOut = NULL;
    if (!shapePath || !shapePath[0])
	return NULL;

    SoGroup *rootGroup = scene_root_group(sceneRoot);
    if (!rootGroup)
	return NULL;
    return scene_shape_find_in_group(rootGroup, shapePath, parentOut);
}

static SbString
scene_shape_owner_source_path(const SoNode *node)
{
    if (!node)
	return "";
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<const SoBRLVListShape *>(node)->
	       ownerSourcePath.getValue();
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<const SoBRLMeshShape *>(node)->
	       ownerSourcePath.getValue();
    return "";
}

static SbString
scene_shape_owner_source_instance_key(const SoNode *node)
{
    if (!node)
	return "";
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<const SoBRLVListShape *>(node)->
	       ownerSourceInstanceKey.getValue();
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<const SoBRLMeshShape *>(node)->
	       ownerSourceInstanceKey.getValue();
    return "";
}

BObolSceneGroupPublishState::BObolSceneGroupPublishState(void) :
    intentPath(nullptr), drawMode(BOBOL_LOD_DRAW_UNKNOWN), fallbackDrawMode(BOBOL_LOD_DRAW_WIRE),
    overlayIntent(FALSE), revalidationRevision(0), visible(TRUE), selected(FALSE), highlighted(FALSE),
    lineStyle(0), lineWidth(0), transparency(0.0f), colorOverride(FALSE), color(1.0f, 1.0f, 1.0f),
    materialColorValid(FALSE), materialColor(1.0f, 1.0f, 1.0f), materialRevision(0)
{}

BObolSceneSummary::BObolSceneSummary(void) :
    valid(FALSE),
    hasRoot(FALSE),
    rootIsGroup(FALSE),
    rootChildCount(0),
    databaseSourceCount(0),
    nonDatabaseRootChildCount(0),
    structuralRevision(0),
    frameRevision(0),
    lastVisitedSourceCount(0),
    lastRealizedSourceCount(0),
    lastFailedSourceCount(0),
    lastDiagnostics("")
{
}

using SceneChildOrders = std::map<SoGroup *, std::vector<SoNode *>>;

/* Only changed parent relations are indexed for prospective reachability. The
 * unchanged graph continues to use native parent auditors. */
struct SceneGraphChanges {
    explicit SceneGraphChanges(const SceneChildOrders &children)
    {
	for (const auto &entry : children) {
	    auto *parent = entry.first;
	    bool changed = size_t(parent->getNumChildren()) != entry.second.size();
	    if (!changed) for (int i = 0; i < parent->getNumChildren(); ++i)
		if (parent->getChild(i) != entry.second[size_t(i)]) { changed = true; break; }
	    if (!changed) continue;
	    this->parents.push_back(parent);
	    std::unordered_set<SoNode *> previous, next(entry.second.begin(), entry.second.end());
	    for (int i = 0; i < parent->getNumChildren(); ++i) previous.insert(parent->getChild(i));
	    for (auto *node : previous) if (!next.count(node)) this->removed[node].insert(parent);
	    for (auto *node : next) if (!previous.count(node)) this->added[node].insert(parent);
	}
    }
    std::vector<SoGroup *> parents;
    std::unordered_map<SoNode *, std::set<SoGroup *>> added, removed;
};

class BObolSceneController::RootOwnership : public SoNodeSensor {
public:
    RootOwnership(BObolSceneController &owner, SoNode *root) : SoNodeSensor(observed, nullptr), scene(owner)
    {
	this->setPriority(0);
	if (root) this->attach(root);
    }
    static std::set<BObolSceneController *> controllers(SoNode *node, const SceneGraphChanges *changes = nullptr);
    template <typename Visit> static void visitControllers(SoNode *node, const SceneGraphChanges *changes, Visit visit);
private:
    // The native auditor identifies a controller at its current root. Scene
    // publications commit peer effects before notification; observing the root
    // does not create another revision effect or schedule policy work.
    static void observed(void *, SoSensor *) {}
    BObolSceneController &scene;
};

struct SceneIndexes {
    std::unordered_map<std::string, SoGroup *> groupPathIndex;
    std::unordered_map<std::string, SoBRLDatabaseSource *> databaseSourcePathIndex;
    std::unordered_map<std::string, std::vector<SoBRLDatabaseSource *>> databaseSourcePathInstancesIndex;
    std::unordered_map<std::string, SoBRLDatabaseSource *> databaseSourceInstanceIndex;
    std::unordered_map<std::string, SoGroup *> databaseSourceInstanceParentIndex;
    std::unordered_map<uint64_t, SoBRLDatabaseSource *> databaseSourceRoutingIndex;
    std::vector<SoBRLDatabaseSource *> databaseSourceOrder;
    std::unordered_map<SoBRLDatabaseSource *, size_t> databaseSourceOrderIndex;

    explicit SceneIndexes(SoNode *root = nullptr)
    {
	if (!root || !root->isOfType(SoGroup::getClassTypeId())) return;
	index_database_sources_recursive(static_cast<SoGroup *>(root), this->groupPathIndex, this->databaseSourcePathIndex,
	    this->databaseSourcePathInstancesIndex, this->databaseSourceInstanceIndex, this->databaseSourceInstanceParentIndex,
	    this->databaseSourceOrder, this->databaseSourceOrderIndex);
	for (auto *source : this->databaseSourceOrder) this->databaseSourceRoutingIndex.emplace(source->getCompactSourceRoutingId(), source);
    }
    void swap(SceneIndexes &other) noexcept
    {
	this->groupPathIndex.swap(other.groupPathIndex);
	this->databaseSourcePathIndex.swap(other.databaseSourcePathIndex);
	this->databaseSourcePathInstancesIndex.swap(other.databaseSourcePathInstancesIndex);
	this->databaseSourceInstanceIndex.swap(other.databaseSourceInstanceIndex);
	this->databaseSourceInstanceParentIndex.swap(other.databaseSourceInstanceParentIndex);
	this->databaseSourceRoutingIndex.swap(other.databaseSourceRoutingIndex);
	this->databaseSourceOrder.swap(other.databaseSourceOrder);
	this->databaseSourceOrderIndex.swap(other.databaseSourceOrderIndex);
    }
};

struct BObolSceneController::Impl {
    Impl(void) :
	root(NULL),
	structuralRevision(0),
	frameRevision(0),
	mutationBatchDepth(0),
	mutationBatchStructuralRevisionPending(FALSE),
	mutationBatchFrameRevisionPending(FALSE),
	databaseSourceIndexValid(FALSE),
	lastVisitedSourceCount(0),
	lastRealizedSourceCount(0),
	lastFailedSourceCount(0),
	lastDiagnostics(""),
	realizationRepository(std::make_shared<BObolRealizationRepository>())
    {
    }

    SoNode *root;
    uint64_t structuralRevision;
    uint64_t frameRevision;
    int mutationBatchDepth;
    SbBool mutationBatchStructuralRevisionPending;
    SbBool mutationBatchFrameRevisionPending;
    mutable SbBool databaseSourceIndexValid;
    SceneIndexes indexes;
    std::unique_ptr<RootOwnership> rootOwnership;
    bool closing = false;
    unsigned int lastVisitedSourceCount;
    unsigned int lastRealizedSourceCount;
    unsigned int lastFailedSourceCount;
    SbString lastDiagnostics;
    SoBRLRealizeAction *activeRealizationAction = nullptr;
    std::shared_ptr<BObolRealizationRepository> realizationRepository;
};

template <typename Visit>
void
BObolSceneController::RootOwnership::visitControllers(SoNode *node, const SceneGraphChanges *changes, Visit visit)
{
    SbList<SoNode *> pending;
    pending.append(node);
    SceneNodeSet visited;
    while (pending.getLength() > 0) {
	auto *current = pending.pop();
	if (!current || !visited.insert(current)) continue;
	if (changes) {
	    auto added = changes->added.find(current);
	    if (added != changes->added.end()) for (auto *parent : added->second) pending.append(parent);
	}
	const auto &auditors = current->getAuditors();
	for (int i = 0; i < auditors.getLength(); ++i) {
	    if (auditors.getType(i) == SoNotRec::PARENT) {
		auto *parent = static_cast<SoGroup *>(auditors.getObject(i));
		if (changes) {
		    auto removed = changes->removed.find(current);
		    if (removed != changes->removed.end() && removed->second.count(parent)) continue;
		}
		pending.append(parent);
	    } else if (auditors.getType(i) == SoNotRec::SENSOR) {
		auto *owner = dynamic_cast<RootOwnership *>(static_cast<SoDataSensor *>(auditors.getObject(i)));
		if (owner && owner->scene.d->root == current && !owner->scene.d->closing) visit(&owner->scene);
	    }
	}
    }
}

std::set<BObolSceneController *>
BObolSceneController::RootOwnership::controllers(SoNode *node, const SceneGraphChanges *changes)
{
    std::set<BObolSceneController *> result;
    visitControllers(node, changes, [&result](BObolSceneController *owner) { result.insert(owner); });
    return result;
}

/* Compose peer revisions across all changed nodes before callbacks. A structural
 * effect includes a frame effect, so an enclosing publication advances each
 * affected controller once even when fields and hierarchy both change. */
class BObolSceneController::PeerEffects {
public:
    enum class Revision { Frame, Structural };
    explicit PeerEffects(BObolSceneController &owner) : scene(owner) {}
    PeerEffects(BObolSceneController &owner, SoNode *target) : PeerEffects(owner) { this->include(target, Revision::Frame); }
    void include(SoNode *target, Revision revision)
    {
	RootOwnership::visitControllers(target, nullptr, [this, revision](BObolSceneController *controller) {
	    if (controller == &this->scene) return;
	    auto inserted = this->peers.emplace(controller, revision);
	    if (revision == Revision::Structural) inserted.first->second = revision;
	});
    }
    void commit() noexcept
    {
	for (const auto &entry : this->peers) {
	    if (entry.second == Revision::Structural) {
		entry.first->clearDatabaseSourceIndex();
		entry.first->advanceStructuralRevision();
	    } else entry.first->advanceFrameRevision();
	}
    }
    static void frameCommitted(void *context) noexcept
    {
	auto &self = *static_cast<PeerEffects *>(context);
	self.commit();
	self.scene.advanceFrameRevision();
    }
private:
    BObolSceneController &scene;
    std::map<BObolSceneController *, Revision> peers;
};

static std::vector<SoBRLDatabaseSource *>
scene_sources(SoNode *root)
{
    std::vector<SoBRLDatabaseSource *> sources;
    std::unordered_set<SoBRLDatabaseSource *> visited;
    if (!scene_walk(root, nullptr, [&](SoNode *node, SoGroup *parent) {
	if (!parent || !node->isOfType(SoBRLDatabaseSource::getClassTypeId())) return;
	auto *source = static_cast<SoBRLDatabaseSource *>(node);
	if (visited.insert(source).second) sources.push_back(source);
    }, [](SoGroup *) -> const std::vector<SoNode *> * { return nullptr; }))
	throw std::invalid_argument("scene root contains a cycle");
    return sources;
}

class BObolSceneController::RootPublication {
public:
    RootPublication(BObolSceneController &owner, SoNode *root) : scene(owner),
	previous(owner.d->root), next(root), repository(owner.d->realizationRepository), indexes(root)
    {
	if (root) this->ownership = std::make_unique<RootOwnership>(owner, root);
	std::vector<BObolRealizationRepository::SourceMembership> memberships;
	const auto oldSources = scene_sources(owner.d->root);
	memberships.reserve(oldSources.size() + this->indexes.databaseSourceOrder.size());
	for (auto *source : oldSources) memberships.push_back({&owner, source, false});
	for (auto *source : this->indexes.databaseSourceOrder) memberships.push_back({&owner, source, true});
	this->sources = this->repository->prepareSourceMembership(memberships, this->indexes.databaseSourceOrder);
    }
    void commit() noexcept
    {
	auto &state = *this->scene.d;
	if (this->next.get()) this->next.get()->ref();
	auto *old = std::exchange(state.root, this->next.get());
	state.rootOwnership.swap(this->ownership);
	state.indexes.swap(this->indexes);
	state.databaseSourceIndexValid = TRUE;
	this->sources->commit();
	this->scene.advanceStructuralRevision();
	// The prepared reference keeps retirement outside the commit.
	if (old) old->unref();
    }
private:
    BObolSceneController &scene;
    SbModernUtils::SoNodeRef previous, next;
    std::shared_ptr<BObolRealizationRepository> repository;
    SceneIndexes indexes;
    std::unique_ptr<RootOwnership> ownership;
    std::unique_ptr<BObolRealizationRepository::SourceUpdate> sources;
};

BObolSceneController::BObolSceneController(void) : d(new Impl)
{
    this->d->realizationRepository->attachController();
}

BObolSceneController::BObolSceneController(SoNode *sceneRoot) : BObolSceneController()
{
    this->setSceneRoot(sceneRoot);
}

BObolSceneController::~BObolSceneController(void)
{
    this->d->closing = true;
    auto *old = std::exchange(this->d->root, nullptr);
    this->d->rootOwnership.reset();
    this->clearDatabaseSourceIndex();
    if (old) this->advanceStructuralRevision();
    this->d->realizationRepository->detachController(this);
    if (old) old->unref();
}

void
BObolSceneController::setSceneRoot(SoNode *sceneRoot)
{
    // A deletion observer cannot resurrect a controller being destroyed.
    if (this->d->closing || this->d->root == sceneRoot) return;
    RootPublication publication(*this, sceneRoot);
    publication.commit();
}

void
BObolSceneController::shareRealizationRepository(BObolSceneController *source)
{
    if (this->d->closing || !source || source == this || !source->d->realizationRepository ||
	this->d->realizationRepository == source->d->realizationRepository) return;
    auto previous = this->d->realizationRepository;
    auto next = source->d->realizationRepository;
    const auto sources = scene_sources(this->d->root);
    std::vector<BObolRealizationRepository::SourceMembership> memberships;
    memberships.reserve(sources.size());
    for (auto *node : sources) memberships.push_back({this, node, false});
    auto retired = previous->prepareSourceMembership(memberships, {});
    for (auto &membership : memberships) membership.present = true;
    auto attached = next->prepareSourceMembership(memberships, sources);
    if (this->d->activeRealizationAction) this->d->activeRealizationAction->stopSceneTraversal();
    this->d->realizationRepository = next;
    retired->commit(); attached->commit();
    next->attachController();
    previous->detachController(this);
}

void
BObolSceneController::clearRealizationRepository(void)
{
    auto repository = this->d->realizationRepository;
    if (!repository)
	return;
    if (this->d->activeRealizationAction)
	this->d->activeRealizationAction->stopSceneTraversal();
    repository->clear();
}

void
BObolSceneController::invalidateRealizationViewVariants(void)
{
    auto repository = this->d->realizationRepository;
    if (!repository)
	return;
    if (this->d->activeRealizationAction)
	this->d->activeRealizationAction->stopSceneTraversal();
    repository->invalidateViewVariants();
}

void
BObolSceneController::renameRealizationObject(
    const char *oldName, const char *newName)
{
    if (!oldName || !oldName[0] || !newName || !newName[0])
	return;
    const std::string oldObject = bobol_realization_cache_object_name(oldName);
    const std::string newObject = bobol_realization_cache_object_name(newName);
    if (oldObject.empty() || newObject.empty() || oldObject == newObject)
	return;
    auto repository = this->d->realizationRepository;
    if (!repository)
	return;
    if (this->d->activeRealizationAction)
	this->d->activeRealizationAction->stopSceneTraversal();
    repository->renameObject(oldObject.c_str(), newObject.c_str());
}

SoNode *
BObolSceneController::getSceneRoot(void) const
{
    return this->d->root;
}

uint64_t
BObolSceneController::getStructuralRevision(void) const
{
    return this->d->structuralRevision;
}

uint64_t
BObolSceneController::getFrameRevision(void) const
{
    return this->d->frameRevision;
}

SbBool
BObolSceneController::getSceneSummary(BObolSceneSummary &summary) const
{
    summary = BObolSceneSummary();
    summary.valid = TRUE;
    summary.hasRoot = this->d->root ? TRUE : FALSE;
    summary.structuralRevision = this->d->structuralRevision;
    summary.frameRevision = this->d->frameRevision;
    summary.lastVisitedSourceCount = this->getLastVisitedSourceCount();
    summary.lastRealizedSourceCount = this->getLastRealizedSourceCount();
    summary.lastFailedSourceCount = this->getLastFailedSourceCount();
    summary.lastDiagnostics = this->getLastDiagnostics();

    if (!this->d->root)
	return TRUE;

    summary.rootIsGroup =
	this->d->root->isOfType(SoGroup::getClassTypeId()) ? TRUE : FALSE;
    if (!summary.rootIsGroup)
	return TRUE;

    SoGroup *group = static_cast<SoGroup *>(this->d->root);
    summary.rootChildCount = group->getNumChildren();
    summary.databaseSourceCount = this->getDatabaseSourceCount();
    int rootDatabaseSourceCount = 0;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *child = group->getChild(i);
	if (child && child->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	    rootDatabaseSourceCount++;
    }
    summary.nonDatabaseRootChildCount =
	summary.rootChildCount - rootDatabaseSourceCount;
    return TRUE;
}

void
BObolSceneController::clearDatabaseSourceIndex(void) const
{
    SceneIndexes empty;
    this->d->indexes.swap(empty);
    this->d->databaseSourceIndexValid = FALSE;
}

void
BObolSceneController::rebuildDatabaseSourceIndex(void) const
{
    BObolPerformanceTimer timer(BOBOL_PERF_SOURCE_INDEX_REBUILD_US);
    if (timer.active()) bobol_performance_counter_add(BOBOL_PERF_SOURCE_INDEX_REBUILD_CALLS, 1);
    SceneIndexes next(this->d->root);
    this->d->indexes.swap(next);
    this->d->databaseSourceIndexValid = TRUE;
}

SoGroup *
BObolSceneController::findIndexedGroup(const char *groupPath) const
{
    if (!groupPath)
	return NULL;

    std::vector<std::string> components;
    scene_group_path_components(groupPath, components);
    if (components.empty())
	return scene_root_group(this->d->root);

    if (!this->d->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();

    const char *normalized = skip_leading_slash(groupPath);
    auto it = this->d->indexes.groupPathIndex.find(
		  std::string(normalized ? normalized : groupPath));
    if (it != this->d->indexes.groupPathIndex.end())
	return it->second;
    return NULL;
}

SoBRLDatabaseSource *
BObolSceneController::findIndexedDatabaseSource(const char *sourcePath) const
{
    if (!sourcePath || !sourcePath[0])
	return NULL;
    if (!this->d->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();

    const char *normalized = skip_leading_slash(sourcePath);
    auto it = this->d->indexes.databaseSourcePathIndex.find(
		  std::string(normalized ? normalized : sourcePath));
    if (it != this->d->indexes.databaseSourcePathIndex.end())
	return it->second;
    return NULL;
}

SbString
BObolSceneController::databaseSourceInstanceKeyForPath(
    const char *sourcePath) const
{
    SoBRLDatabaseSource *source = this->findIndexedDatabaseSource(sourcePath);
    return database_source_effective_instance_key(source);
}

SoBRLDatabaseSource *
BObolSceneController::findIndexedDatabaseSourceInstance(
    const char *sourceInstanceKey) const
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return NULL;
    if (!this->d->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();

    std::string normalizedKey = normalized_index_key(sourceInstanceKey);
    auto it = this->d->indexes.databaseSourceInstanceIndex.find(normalizedKey);
    if (it != this->d->indexes.databaseSourceInstanceIndex.end())
	return it->second;
    return NULL;
}

SoGroup *
BObolSceneController::findIndexedDatabaseSourceInstanceParent(
    const char *sourceInstanceKey) const
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return NULL;
    if (!this->d->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();

    std::string normalizedKey = normalized_index_key(sourceInstanceKey);
    auto it = this->d->indexes.databaseSourceInstanceParentIndex.find(normalizedKey);
    if (it != this->d->indexes.databaseSourceInstanceParentIndex.end())
	return it->second;
    return NULL;
}

static int
scene_retire_compact_hierarchy_descendants(BObolSceneController *scene,
	const char *compactRootInstanceKey)
{
    if (!scene)
	return 0;

    struct Source {
	SbModernUtils::SoNodeRef owner;
	SbString key;
	SbString parent;
	std::vector<size_t> ancestors;
	bool compact = false;
	bool descendant = false;
    };
    std::vector<Source> sources;
    const int sourceCount = scene->getDatabaseSourceCount();
    sources.reserve(size_t(sourceCount));
    for (int i = 0; i < sourceCount; i++) {
	SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	BObolDatabaseSourceSummary summary;
	if (!source || !source->getSummary(summary) || !summary.valid ||
	    summary.instanceKey.getLength() == 0) continue;
	sources.push_back({SbModernUtils::SoNodeRef(source), summary.instanceKey,
	    summary.parentInstanceKey, {}, bool(source->hasCompactInstanceIndex()), false});
    }
    const size_t ambiguous = std::numeric_limits<size_t>::max();
    std::unordered_map<std::string, size_t> identities;
    for (size_t i = 0; i < sources.size(); ++i) {
	auto inserted = identities.emplace(sources[i].key.getString(), i);
	if (!inserted.second) inserted.first->second = ambiguous;
    }
    size_t compactRoot = ambiguous;
    if (compactRootInstanceKey && compactRootInstanceKey[0]) {
	const auto root = identities.find(compactRootInstanceKey);
	if (root == identities.end() || root->second == ambiguous)
	    return 0;
	compactRoot = root->second;
    }

    bool found = true;
    while (found) {
	found = false;
	for (auto &source : sources) {
	    if (source.compact || source.descendant || source.parent.getLength() == 0) continue;
	    const auto parent = identities.find(source.parent.getString());
	    if (parent == identities.end() || parent->second == ambiguous) continue;
	    const auto &owner = sources[parent->second];
	    if (!owner.compact && !owner.descendant) continue;
	    source.ancestors = owner.ancestors;
	    source.ancestors.push_back(parent->second);
	    source.descendant = true;
	    found = true;
	}
    }

    const auto current = [scene](const Source &source) {
	auto *node = static_cast<SoBRLDatabaseSource *>(source.owner.get());
	return node && scene->findDatabaseSourceInstance(source.key.getString()) == node &&
	    database_source_effective_instance_key(node) == source.key &&
	    node->parentInstanceKey.getValue() == source.parent;
    };
    int retired = 0;
    for (auto target = sources.rbegin(); target != sources.rend(); ++target) {
	if (!target->descendant || !current(*target)) continue;
	if (compactRoot != ambiguous &&
	    std::find(target->ancestors.begin(), target->ancestors.end(),
		compactRoot) == target->ancestors.end())
	    continue;
	bool accepted = true;
	for (size_t ancestor : target->ancestors) {
	    const auto &owner = sources[ancestor];
	    accepted = accepted && current(owner) &&
		(!owner.compact || static_cast<SoBRLDatabaseSource *>(owner.owner.get())->hasCompactInstanceIndex());
	}
	if (accepted && scene->removeDatabaseSourceInstance(target->key.getString()) > 0) retired = 1;
    }
    return retired;
}

SbBool
BObolSceneController::realizePending(void)
{
    return this->realizePending(nullptr);
}

SbBool
BObolSceneController::realizePending(BObolSourceRealizationEffects *effects)
{
    return this->realizeSubtree(this->d->root, effects, nullptr);
}

SbBool
BObolSceneController::realizeDatabaseSourceInstance(
    const char *sourceInstanceKey,
    const BObolSourceRealizationStamp &stamp)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return FALSE;
    SoBRLDatabaseSource *source =
	this->findIndexedDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return FALSE;
    SbModernUtils::SoNodeRef sourceOwner(source);
    if (this->findIndexedDatabaseSourceInstance(sourceInstanceKey) != source ||
	!source->matchesRealizationStamp(stamp))
	return FALSE;
    return this->realizeSubtree(source, nullptr, sourceInstanceKey);
}

SbBool
BObolSceneController::adoptDatabaseSourceInstanceMeshLod(
    const char *sourceInstanceKey,
    const BObolSourceRealizationStamp &stamp,
    struct BObolMeshLod *lod,
    const SbVec3f &bmin,
    const SbVec3f &bmax)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0] || !lod)
	return FALSE;
    SoBRLDatabaseSource *source =
	this->findIndexedDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return FALSE;
    SbModernUtils::SoNodeRef sourceOwner(source);
    if (this->findIndexedDatabaseSourceInstance(sourceInstanceKey) != source ||
	!source->matchesRealizationStamp(stamp))
	return FALSE;
    return source->adoptMeshLod(stamp, lod, bmin, bmax) > 0 ? TRUE : FALSE;
}

SbBool
BObolSceneController::realizeSubtree(SoNode *root,
    BObolSourceRealizationEffects *effects,
    const char *compactRootInstanceKey)
{
    this->d->lastVisitedSourceCount = 0;
    this->d->lastRealizedSourceCount = 0;
    this->d->lastFailedSourceCount = 0;
    this->d->lastDiagnostics = "";

    if (!root)
	return FALSE;

    /* Callbacks may replace either owner. Keep the captured inputs alive until
     * the action and its completed progress have finished unwinding. */
    const auto repository = this->d->realizationRepository;
    SbModernUtils::SoNodeRef rootOwner(root);
    auto publication = prepareRealizationEffects(*this, effects);
    SoBRLRealizeAction action;
    action.setRealizationRepository(repository.get());
    action.setRetainRealizationCache(TRUE);
    action.publicationEffects = publication.get();

    /* Observers read the action's committed progress without a second string
     * allocation. Preserve it on every exit, including a later traversal or
     * repository failure; nested calls restore the enclosing action. */
    struct ProgressScope {
	ProgressScope(Impl &owner, SoBRLRealizeAction &current) :
	    state(owner), action(current)
	{
	    action.enclosingSceneAction = state.activeRealizationAction;
	    state.activeRealizationAction = &current;
	}
	~ProgressScope()
	{
	    state.lastVisitedSourceCount = action.getVisitedSourceCount();
	    state.lastRealizedSourceCount = action.getRealizedSourceCount();
	    state.lastFailedSourceCount = action.getFailedSourceCount();
	    state.lastDiagnostics = std::move(action.diagnostics);
	    state.activeRealizationAction = action.enclosingSceneAction;
	    action.enclosingSceneAction = nullptr;
	}
	Impl &state;
	SoBRLRealizeAction &action;
    } progress(*this->d, action);
    action.apply(rootOwner.get());
    if (!action.hasTerminated())
	(void)scene_retire_compact_hierarchy_descendants(this,
	    compactRootInstanceKey);
    return action.getFailedSourceCount() == 0;
}

void
BObolSceneController::beginSceneMutationBatch(size_t expectedDatabaseSources,
	size_t expectedGroups)
{
    if (this->d->mutationBatchDepth == 0) {
	this->d->mutationBatchStructuralRevisionPending = FALSE;
	this->d->mutationBatchFrameRevisionPending = FALSE;
	if (this->d->root && this->d->root->isOfType(SoGroup::getClassTypeId()) &&
	    !this->d->databaseSourceIndexValid)
	    this->rebuildDatabaseSourceIndex();
	if (expectedDatabaseSources > 0) {
	    this->d->indexes.databaseSourcePathIndex.reserve(
		this->d->indexes.databaseSourcePathIndex.size() + expectedDatabaseSources);
	    this->d->indexes.databaseSourcePathInstancesIndex.reserve(
		this->d->indexes.databaseSourcePathInstancesIndex.size() +
		expectedDatabaseSources);
	    this->d->indexes.databaseSourceInstanceIndex.reserve(
		this->d->indexes.databaseSourceInstanceIndex.size() +
		expectedDatabaseSources);
	    this->d->indexes.databaseSourceInstanceParentIndex.reserve(
		this->d->indexes.databaseSourceInstanceParentIndex.size() +
		expectedDatabaseSources);
	    this->d->indexes.databaseSourceOrder.reserve(
		this->d->indexes.databaseSourceOrder.size() + expectedDatabaseSources);
	    this->d->indexes.databaseSourceOrderIndex.reserve(
		this->d->indexes.databaseSourceOrderIndex.size() + expectedDatabaseSources);
	}
	if (expectedGroups > 0)
	    this->d->indexes.groupPathIndex.reserve(
		this->d->indexes.groupPathIndex.size() + expectedGroups);
    }

    this->d->mutationBatchDepth++;
}

void
BObolSceneController::endSceneMutationBatch(void)
{
    if (this->d->mutationBatchDepth <= 0)
	return;

    this->d->mutationBatchDepth--;
    if (this->d->mutationBatchDepth > 0)
	return;

    const SbBool structuralChanged =
	this->d->mutationBatchStructuralRevisionPending;
    const SbBool frameChanged =
	this->d->mutationBatchFrameRevisionPending;
    this->d->mutationBatchStructuralRevisionPending = FALSE;
    this->d->mutationBatchFrameRevisionPending = FALSE;

    if (structuralChanged) {
	bobol_identity_advance(this->d->structuralRevision);
	bobol_identity_advance(this->d->frameRevision);
    } else if (frameChanged) {
	bobol_identity_advance(this->d->frameRevision);
    }
}

SoGroup *
BObolSceneController::findGroup(const char *groupPath) const
{
    return this->findIndexedGroup(groupPath);
}

/* Reserve hash capacity and detached nodes while the old entries remain live.
 * Commit changes only the selected keys and never allocates or notifies. */
template <typename Value>
class PreparedSceneIndex {
public:
    using Map = std::unordered_map<std::string, Value *>;
    explicit PreparedSceneIndex(Map &index) : target(index) {}
    void set(const std::string &key, Value *value)
    {
	if (!key.empty()) this->changes[key] = value;
    }
    bool contains(const std::string &key) const { return this->changes.count(key) != 0; }
    bool empty() const { return this->changes.empty(); }
    void prepare()
    {
	size_t added = 0;
	for (const auto &entry : this->changes)
	    if (entry.second && this->target.find(entry.first) == this->target.end()) {
		this->nodes.emplace(entry);
		++added;
	    }
	if (added) this->target.reserve(this->target.size() + added);
    }
    void commit() noexcept
    {
	for (const auto &entry : this->changes) {
	    auto found = this->target.find(entry.first);
	    if (!entry.second) {
		if (found != this->target.end()) this->target.erase(found);
	    } else if (found != this->target.end()) found->second = entry.second;
	    else this->target.insert(this->nodes.extract(entry.first));
	}
    }
private:
    Map &target;
    Map changes, nodes;
};

class BObolSceneController::MembershipPublication {
public:
    MembershipPublication(PeerEffects &effects, const SceneChildOrders &children,
	PeerEffects::Revision revision)
    {
	SceneGraphChanges graph(children);
	std::vector<SoNode *> pending;
	for (const auto &entry : graph.added) pending.push_back(entry.first);
	for (const auto &entry : graph.removed) pending.push_back(entry.first);
	std::unordered_set<SoNode *> visited;
	std::unordered_set<SoBRLDatabaseSource *> sources;
	while (!pending.empty()) {
	    auto *node = pending.back(); pending.pop_back();
	    if (!node || !visited.insert(node).second) continue;
	    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) sources.insert(static_cast<SoBRLDatabaseSource *>(node));
	    if (!node->isOfType(SoGroup::getClassTypeId())) continue;
	    auto *group = static_cast<SoGroup *>(node);
	    for (int i = 0; i < group->getNumChildren(); ++i) pending.push_back(group->getChild(i));
	    auto next = children.find(group);
	    if (next != children.end()) pending.insert(pending.end(), next->second.begin(), next->second.end());
	}
	this->sourceTopology = !sources.empty();
	const PeerEffects::Revision effectiveRevision = this->sourceTopology ?
	    PeerEffects::Revision::Structural : revision;
	for (auto *parent : graph.parents)
	    effects.include(parent, effectiveRevision);
	struct Draft {
	    std::shared_ptr<BObolRealizationRepository> repository;
	    std::vector<BObolRealizationRepository::SourceMembership> memberships;
	    std::set<SoBRLDatabaseSource *> seeded;
	};
	std::map<BObolRealizationRepository *, Draft> drafts;
	auto stage = [&drafts](BObolSceneController *controller, SoBRLDatabaseSource *source, bool present) {
	    auto repository = controller->d->realizationRepository;
	    if (repository->hasSourceOwner(controller, source) == present) return;
	    auto &draft = drafts[repository.get()];
	    draft.repository = repository;
	    draft.memberships.push_back({controller, source, present});
	    if (present) draft.seeded.insert(source);
	};
	for (auto *source : sources) {
	    const auto previous = RootOwnership::controllers(source);
	    const auto next = RootOwnership::controllers(source, &graph);
	    for (auto *controller : previous)
		if (controller->d->root != source && !next.count(controller)) stage(controller, source, false);
	    for (auto *controller : next)
		if (controller->d->root != source) stage(controller, source, true);
	}
	for (auto &entry : drafts) {
	    auto &draft = entry.second;
	    auto update = draft.repository->prepareSourceMembership(draft.memberships,
		std::vector<SoBRLDatabaseSource *>(draft.seeded.begin(), draft.seeded.end()));
	    this->updates.push_back({std::move(draft.repository), std::move(update)});
	}
    }
    void commit() noexcept
    {
	for (auto &update : this->updates) update.fields->commit();
    }
    bool changesSourceTopology() const noexcept { return this->sourceTopology; }
private:
    struct Update {
	std::shared_ptr<BObolRealizationRepository> repository;
	std::unique_ptr<BObolRealizationRepository::SourceUpdate> fields;
    };
    std::vector<Update> updates;
    bool sourceTopology = false;
};

class BObolSceneController::HierarchyPublication {
public:
    explicit HierarchyPublication(BObolSceneController &owner) : scene(owner), peers(owner), groups(owner.d->indexes.groupPathIndex)
    {
	this->retainedNodes.emplace_back(owner.d->root);
	if (!owner.d->databaseSourceIndexValid) owner.rebuildDatabaseSourceIndex();
    }
    SoGroup *ensurePath(const char *path)
    {
	SoGroup *target = scene_root_group(this->scene.d->root);
	std::vector<std::string> components;
	scene_group_path_components(path, components);
	for (const auto &component : components) {
	    this->groupPath = scene_group_append_path(this->groupPath, component.c_str());
	    auto *group = this->scene.findIndexedGroup(this->groupPath.getString());
	    if (!group) {
		auto *created = new SoBRLSceneGroup;
		SbModernUtils::SoNodeRef retained(created);
		this->retainedNodes.push_back(std::move(retained));
		created->setName(SbName(component.c_str()));
		created->groupPath = this->groupPath;
		this->nextChildren(target).push_back(created);
		this->groups.set(normalized_index_key(this->groupPath.getString()), created);
		group = created;
	    }
	    target = group;
	}
	return target;
    }
    bool reaches(SoGroup *ancestor, SoGroup *target) const
    {
	std::vector<SoGroup *> pending{ancestor};
	std::unordered_set<SoGroup *> visited;
	while (!pending.empty()) {
	    auto *group = pending.back(); pending.pop_back();
	    if (group == target) return true;
	    if (!visited.insert(group).second) continue;
	    auto append = [&pending](SoNode *child) {
		if (child->isOfType(SoGroup::getClassTypeId())) pending.push_back(static_cast<SoGroup *>(child));
	    };
	    const auto prepared = this->childOrders.find(group);
	    if (prepared != this->childOrders.end()) {
		for (auto *child : prepared->second) append(child);
	    } else {
		for (int i = 0; i < group->getNumChildren(); ++i) append(group->getChild(i));
	    }
	}
	return false;
    }
    const std::vector<SoNode *> *preparedChildren(SoGroup *group) const
    {
	auto found = this->childOrders.find(group);
	return found == this->childOrders.end() ? nullptr : &found->second;
    }
    PeerEffects &effects() { return this->peers; }
    const SbString &path() const { return this->groupPath; }
    bool changed() const { return !this->childOrders.empty(); }
    bool changesSourceTopology() const noexcept
    {
	return this->membership && this->membership->changesSourceTopology();
    }
    std::vector<SoNode *> &nextChildren(SoGroup *parent)
    {
	auto inserted = this->childOrders.emplace(parent, std::vector<SoNode *>());
	auto &order = inserted.first->second;
	if (inserted.second) {
	    this->retainedNodes.emplace_back(parent);
	    order.reserve(size_t(parent->getNumChildren()) + 1);
	    for (int i = 0; i < parent->getNumChildren(); ++i) order.push_back(parent->getChild(i));
	}
	return order;
    }
    void removeChild(SoGroup *parent, int index)
    {
	auto &order = this->nextChildren(parent);
	this->removals.emplace(parent, parent->getChildren()->prepareRemoval({index}));
	order.erase(order.begin() + index);
    }
    void prepare(bool prepareChildren = true,
	PeerEffects::Revision revision = PeerEffects::Revision::Structural)
    {
	if (prepareChildren) {
	    for (auto &edit : this->childOrders)
		if (!this->removals.count(edit.first))
		    this->children.push_back(edit.first->getChildren()->prepareReplacement(edit.second));
	}
	this->groups.prepare();
	if (!this->childOrders.empty())
	    this->membership = std::make_unique<MembershipPublication>(
		this->peers, this->childOrders, revision);
    }
    void commit() noexcept
    {
	for (auto &edit : this->children) edit->commit();
	for (auto &edit : this->removals) edit.second->commit();
	this->commitEffects();
    }
    void commitEffects() noexcept
    {
	this->groups.commit();
	if (this->membership) this->membership->commit();
	this->peers.commit();
    }
    void notify(std::exception_ptr &failure)
    {
	auto notify = [&failure](auto &edit) {
	    try { edit->notify(); }
	    catch (...) { if (!failure) failure = std::current_exception(); }
	};
	for (auto &edit : this->children) notify(edit);
	for (auto &edit : this->removals) notify(edit.second);
    }
private:
    BObolSceneController &scene;
    PeerEffects peers;
    SbString groupPath;
    std::vector<SbModernUtils::SoNodeRef> retainedNodes;
    SceneChildOrders childOrders;
    std::unique_ptr<MembershipPublication> membership;
    std::vector<std::unique_ptr<SoChildList::Replacement>> children;
    std::map<SoGroup *, std::unique_ptr<SoChildList::Removal>> removals;
    PreparedSceneIndex<SoGroup> groups;
};

SoGroup *
BObolSceneController::ensureGroup(const char *groupPath)
{
    if (!groupPath || !scene_root_group(this->d->root)) return nullptr;
    HierarchyPublication publication(*this);
    SoGroup *target = publication.ensurePath(groupPath);
    if (!publication.changed()) return target;
    publication.prepare();
    publication.commit();
    this->advanceStructuralRevision();
    std::exception_ptr failure;
    publication.notify(failure);
    if (failure) std::rethrow_exception(failure);
    // Observers can replace, rename or remove the prepared target. Resolve
    // the canonical path again instead of returning a detached borrowed node.
    return this->findIndexedGroup(publication.path().getString());
}

class BObolSceneController::GroupPublication {
public:
    GroupPublication(BObolSceneController &owner, SoBRLSceneGroup &target) :
	scene(owner), retained(&target), candidate(new SoBRLSceneGroup)
    {
	this->next().enableNotify(FALSE);
	copy_publication_scalar_fields(this->next(), target);
    }
    void setIntent(const char *groupPath, const char *intentPath, int drawMode,
	int fallbackDrawMode, SbBool overlayIntent, uint32_t revalidationRevision)
    {
	auto &value = this->next();
	SbString retainedPath = value.groupPath.getValue();
	if (retainedPath.getLength() == 0) retainedPath = skip_leading_slash(groupPath);
	value.drawIntentValid = TRUE;
	value.drawIntentPath = intentPath && intentPath[0] ? SbString(intentPath) : retainedPath;
	value.drawMode = drawMode;
	value.fallbackDrawMode = fallbackDrawMode;
	value.overlayIntent = overlayIntent;
	value.revalidationRevision = revalidationRevision;
    }
    void setDisplay(SbBool visible, SbBool selected, SbBool highlighted, int lineStyle,
	int lineWidth, float transparency, SbBool colorOverride, const SbColor &color,
	SbBool materialColorValid, const SbColor &materialColor, uint32_t materialRevision)
    {
	auto &value = this->next();
	value.visible = visible;
	value.selected = selected;
	value.highlighted = highlighted;
	value.lineStyle = lineStyle;
	value.lineWidth = lineWidth;
	if (scene_group_float_different(value.transparency.getValue(), transparency))
	    value.transparency = transparency;
	value.colorOverride = colorOverride;
	if (!scene_group_color_equal(value.color.getValue(), color)) value.color = color;
	value.materialColorValid = materialColorValid;
	if (!scene_group_color_equal(value.materialColor.getValue(), materialColor)) value.materialColor = materialColor;
	value.materialRevision = materialRevision;
    }
    void setDisplayPatch(const BObolDatabaseSourceDisplayPatch &patch)
    {
	auto &value = this->next();
	if (patch.visibleValid) value.visible = patch.visible;
	if (patch.selectedValid) value.selected = patch.selected;
	if (patch.highlightedValid) value.highlighted = patch.highlighted;
	if (patch.lineStyleValid) value.lineStyle = patch.lineStyle;
	if (patch.lineWidthValid) value.lineWidth = patch.lineWidth;
	if (patch.transparencyValid &&
	    scene_group_float_different(value.transparency.getValue(),
		patch.transparency))
	    value.transparency = patch.transparency;
	if (patch.colorOverrideValid)
	    value.colorOverride = patch.colorOverride;
	if (patch.colorValid &&
	    !scene_group_color_equal(value.color.getValue(), patch.color))
	    value.color = patch.color;
    }
    void setPath(const SbString &path) { this->next().groupPath = path; }
    void setIntentPath(const SbString &path)
    {
	if (this->next().drawIntentValid.getValue())
	    this->next().drawIntentPath = path;
    }
    void prepare()
    {
	this->fields = std::make_unique<PreparedScalarFields>(*this->retained.get(), this->next());
    }
    bool changed() const { return this->fields->changed(); }
    void commit() { this->fields->commit(); }
    void restore() { this->fields->restore(); }
    void notify(std::exception_ptr &failure)
    {
	if (this->changed()) this->fields->notify(failure);
    }
    int publish()
    {
	this->prepare();
	if (!this->changed()) return 0;
	PeerEffects effects(this->scene, this->retained.get());
	this->commit();
	effects.commit();
	this->scene.advanceFrameRevision();
	this->restore();
	std::exception_ptr failure;
	this->notify(failure);
	if (failure) std::rethrow_exception(failure);
	return 1;
    }
private:
    SoBRLSceneGroup &next() { return *static_cast<SoBRLSceneGroup *>(this->candidate.get()); }
    BObolSceneController &scene;
    SbModernUtils::SoNodeRef retained, candidate;
    std::unique_ptr<PreparedScalarFields> fields;
};

int
BObolSceneController::setGroupDrawIntent(const char *groupPath, const char *intentPath,
    int drawMode, int fallbackDrawMode, SbBool overlayIntent, uint32_t revalidationRevision)
{
    SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId())) return -1;
    GroupPublication publication(*this, *static_cast<SoBRLSceneGroup *>(group));
    publication.setIntent(groupPath, intentPath, drawMode, fallbackDrawMode, overlayIntent, revalidationRevision);
    return publication.publish();
}

int
BObolSceneController::setGroupDisplayState(const char *groupPath, SbBool visible,
    SbBool selected, SbBool highlighted, int lineStyle, int lineWidth, float transparency,
    SbBool colorOverride, const SbColor &color, SbBool materialColorValid,
    const SbColor &materialColor, uint32_t materialRevision)
{
    SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId())) return -1;
    GroupPublication publication(*this, *static_cast<SoBRLSceneGroup *>(group));
    publication.setDisplay(visible, selected, highlighted, lineStyle, lineWidth, transparency,
	colorOverride, color, materialColorValid, materialColor, materialRevision);
    return publication.publish();
}

int
BObolSceneController::renameGroup(const char *groupPath,
				  const char *newLeafName)
{
    if (!newLeafName || !newLeafName[0] || strchr(newLeafName, '/'))
	return 0;

    std::vector<std::string> components;
    scene_group_path_components(groupPath, components);
    if (components.empty())
	return 0;

    SoGroup *parent = scene_root_group(this->d->root);
    if (!parent)
	return -1;

    for (size_t i = 0; i + 1 < components.size(); i++) {
	parent = scene_group_find_child(parent, components[i].c_str());
	if (!parent)
	    return 0;
    }

    SoGroup *target = scene_group_find_child(parent,
		      components.back().c_str());
    if (!target)
	return 0;
    if (scene_group_node_name_equal(target, newLeafName))
	return 0;
    if (scene_group_find_child(parent, newLeafName))
	return 0;

    const SbString parentPath =
	scene_group_component_path(components, components.size() - 1);
    const SbString newGroupPath =
	scene_group_append_path(parentPath, newLeafName);
    SbModernUtils::SoNodeRef retainedRoot(this->d->root), retainedTarget(target);
    if (!this->d->databaseSourceIndexValid) this->rebuildDatabaseSourceIndex();
    PreparedSceneIndex<SoGroup> groups(this->d->indexes.groupPathIndex);
    std::vector<std::unique_ptr<GroupPublication>> updates;
    PeerEffects effects(*this);
    effects.include(target, PeerEffects::Revision::Structural);
    struct Path { SoGroup *group; SbString value; };
    std::vector<Path> pending{{target, newGroupPath}};
    std::unordered_set<SoGroup *> visited;
    while (!pending.empty()) {
	Path path = std::move(pending.back()); pending.pop_back();
	if (!visited.insert(path.group).second) continue;
	if (path.group->isOfType(SoBRLSceneGroup::getClassTypeId())) {
	    auto &group = *static_cast<SoBRLSceneGroup *>(path.group);
	    groups.set(normalized_index_key(group.groupPath.getValue().getString()), nullptr);
	    groups.set(normalized_index_key(path.value.getString()), &group);
	    auto update = std::make_unique<GroupPublication>(*this, group);
	    update->setPath(path.value); update->prepare();
	    effects.include(&group, PeerEffects::Revision::Structural);
	    updates.push_back(std::move(update));
	}
	// Visit the last child first so a shared group keeps the final path of
	// the preceding recursive traversal, while preparing each object once.
	for (int i = 0; i < path.group->getNumChildren(); ++i) {
	    auto *child = path.group->getChild(i);
	    if (child->isOfType(SoGroup::getClassTypeId()))
		pending.push_back({static_cast<SoGroup *>(child), scene_child_summary_path(path.value, child)});
	}
    }
    groups.prepare();
    // setName preserves the old registry on allocation failure and sends no
    // notifications. Every remaining scene write has already been prepared.
    target->setName(SbName(newLeafName));
    for (auto &update : updates) update->commit();
    groups.commit();
    effects.commit();
    this->advanceStructuralRevision();
    for (auto &update : updates) update->restore();
    std::exception_ptr failure;
    for (auto &update : updates) update->notify(failure);
    if (failure) std::rethrow_exception(failure);
    return 1;
}

class BObolSceneController::ChildPublication {
public:
    enum class Effect { Frame, Structural };

    ChildPublication(BObolSceneController &owner, SoGroup &parent, SoNode &child, bool append) : ChildPublication(owner)
    {
	this->prepareEdit(parent, child, append);
    }
    ChildPublication(BObolSceneController &owner, SoGroup &parent) : ChildPublication(owner)
    {
	auto &children = this->hierarchy.nextChildren(&parent);
	this->retainedChildren.reserve(children.size());
	for (auto *child : children) this->retainedChildren.emplace_back(child);
	children.clear();
	this->prepare(parent, false);
    }
    ChildPublication(BObolSceneController &owner, SoGroup &parent,
	const std::vector<SoNode *> &children,
	Effect requestedEffect = Effect::Structural) : ChildPublication(owner)
    {
	this->effect = requestedEffect;
	for (int i = 0; i < parent.getNumChildren(); ++i) this->retainChild(parent.getChild(i));
	for (auto *child : children) this->retainChild(child);
	this->hierarchy.nextChildren(&parent) = children;
	this->prepare(parent, false, false);
	if (!this->valid)
	    throw std::invalid_argument("realization child replacement contains a cycle");
    }
    ChildPublication(BObolSceneController &owner, const SceneChildOrders &orders,
	bool prepareChildren = false) : ChildPublication(owner)
    {
	if (orders.empty()) throw std::invalid_argument("source child publication has no child orders");
	for (const auto &entry : orders) {
	    auto *parent = entry.first;
	    if (!parent) throw std::invalid_argument("source child publication has no parent");
	    for (int i = 0; i < parent->getNumChildren(); ++i) this->retainChild(parent->getChild(i));
	    for (auto *child : entry.second) this->retainChild(child);
	    this->hierarchy.nextChildren(parent) = entry.second;
	}
	this->prepare(*orders.begin()->first, false, prepareChildren);
	if (!this->valid) throw std::invalid_argument("source child publication contains a cycle");
    }
    void includeEffect(SoNode *node)
    {
	this->hierarchy.effects().include(node, PeerEffects::Revision::Frame);
    }
    void includeStructuralEffect(SoNode *node)
    {
	this->hierarchy.effects().include(node, PeerEffects::Revision::Structural);
    }
    void commitEffects() noexcept { this->commit(false); }
    static int moveShape(BObolSceneController &owner, SoGroup &parent, SoGroup &destination, SoNode &shape)
    {
	ChildPublication publication(owner);
	publication.retainedChildren.emplace_back(&shape);
	publication.hierarchy.removeChild(&parent, scene_group_find_child_index(&parent, &shape));
	publication.hierarchy.nextChildren(&destination).push_back(&shape);
	publication.prepare(parent, false);
	return publication.publish();
    }
    static int removeSource(BObolSceneController &owner, SoBRLDatabaseSource *source)
    {
	if (!source) return 0;
	ChildPublication publication(owner);
	const SbString key = database_source_effective_instance_key(source);
	SoGroup *parent = owner.findIndexedDatabaseSourceInstanceParent(key.getString());
	if (!parent || scene_group_find_child_index(parent, source) < 0) {
	    // A colliding instance key can index another source's parent. Keep
	    // the node selected by the path lookup as the removal target.
	    parent = nullptr;
	    if (!publication.walk(owner.d->root, nullptr, [&](SoNode *node, SoGroup *nodeParent) {
		if (node == source) parent = nodeParent;
	    })) return -1;
	}
	if (!parent) return 0;
	publication.prepareEdit(*parent, *source, false);
	return publication.publish();
    }
    static int clearSources(BObolSceneController &owner, SoGroup &root)
    {
	ChildPublication publication(owner);
	std::vector<SoGroup *> pending{&root};
	std::unordered_set<SoGroup *> visited;
	auto removable = [](SoNode *node) {
	    return node && node->isOfType(SoBRLDatabaseSource::getClassTypeId()) &&
		!static_cast<SoBRLDatabaseSource *>(node)->auxiliarySource.getValue();
	};
	while (!pending.empty()) {
	    auto *parent = pending.back(); pending.pop_back();
	    if (!visited.insert(parent).second) continue;
	    bool changes = false;
	    for (int i = 0; i < parent->getNumChildren(); ++i) {
		auto *child = parent->getChild(i);
		if (!child) continue;
		if (child->isOfType(SoBRLDatabaseSource::getClassTypeId())) changes = changes || removable(child);
		else if (child->isOfType(SoGroup::getClassTypeId())) pending.push_back(static_cast<SoGroup *>(child));
	    }
	    if (!changes) continue;
	    auto &children = publication.hierarchy.nextChildren(parent);
	    children.erase(std::remove_if(children.begin(), children.end(), [&](SoNode *child) {
		if (!removable(child)) return false;
		publication.retainedChildren.emplace_back(child);
		return true;
	    }), children.end());
	}
	const size_t removed = publication.retainedChildren.size();
	if (!removed) return 0;
	if (removed > size_t(std::numeric_limits<int>::max())) return -1;
	publication.prepare(root, false);
	return publication.publish() < 0 ? -1 : int(removed);
    }
    int publish()
    {
	if (!this->valid) return -1;
	this->commitPrepared();
	std::exception_ptr failure;
	this->notify(failure);
	if (failure) std::rethrow_exception(failure);
	return 1;
    }
    void commitPrepared() noexcept { this->commit(true); }
    void notify(std::exception_ptr &failure) { this->hierarchy.notify(failure); }
private:
    void retainChild(SoNode *node)
    {
	if (node && std::none_of(this->retainedChildren.begin(), this->retainedChildren.end(),
	    [node](const auto &current) { return current.get() == node; }))
	    this->retainedChildren.emplace_back(node);
    }
    void commit(bool children) noexcept
    {
	auto &state = *this->scene.d;
	if (children) this->hierarchy.commit();
	else this->hierarchy.commitEffects();
	for (auto it = this->pathValues.begin(); it != this->pathValues.end();) {
	    auto next = it++;
	    auto previous = state.indexes.databaseSourcePathInstancesIndex.find(next->first);
	    if (next->second.empty()) {
		if (previous != state.indexes.databaseSourcePathInstancesIndex.end()) state.indexes.databaseSourcePathInstancesIndex.erase(previous);
	    } else if (previous != state.indexes.databaseSourcePathInstancesIndex.end()) previous->second.swap(next->second);
	    else state.indexes.databaseSourcePathInstancesIndex.insert(this->pathValues.extract(next));
	}
	this->groups.commit(); this->paths.commit(); this->instances.commit(); this->parents.commit();
	for (auto *source : this->retired) {
	    state.indexes.databaseSourceOrderIndex.erase(source);
	    state.indexes.databaseSourceRoutingIndex.erase(source->getCompactSourceRoutingId());
	}
	auto &order = state.indexes.databaseSourceOrder;
	if (!this->retired.empty()) {
	    order.erase(std::remove_if(order.begin(), order.end(), [this](SoBRLDatabaseSource *source) {
		return this->affected.count(source) && !this->reachable.count(source);
	    }), order.end());
	    for (size_t i = 0; i < order.size(); ++i) state.indexes.databaseSourceOrderIndex.find(order[i])->second = i;
	}
	for (auto *source : this->added) {
	    auto index = this->newOrder.extract(source);
	    index.mapped() = order.size();
	    order.push_back(source);
	    state.indexes.databaseSourceOrderIndex.insert(std::move(index));
	    state.indexes.databaseSourceRoutingIndex.insert(this->newRouting.extract(source->getCompactSourceRoutingId()));
	}
	if (this->effect == Effect::Structural ||
	    this->hierarchy.changesSourceTopology())
	    this->scene.advanceStructuralRevision(children ? TRUE : FALSE);
	else
	    this->scene.advanceFrameRevision();
    }
    explicit ChildPublication(BObolSceneController &owner) :
	scene(owner), hierarchy(owner),
	groups(owner.d->indexes.groupPathIndex), paths(owner.d->indexes.databaseSourcePathIndex),
	instances(owner.d->indexes.databaseSourceInstanceIndex), parents(owner.d->indexes.databaseSourceInstanceParentIndex)
    {}
    void prepareEdit(SoGroup &parent, SoNode &child, bool append)
    {
	this->retainedChildren.emplace_back(&child);
	if (append && child.isOfType(SoGroup::getClassTypeId()) &&
	    this->hierarchy.reaches(static_cast<SoGroup *>(&child), &parent)) return;
	if (append) this->hierarchy.nextChildren(&parent).push_back(&child);
	else this->hierarchy.removeChild(&parent, scene_group_find_child_index(&parent, &child));
	this->prepare(parent, append);
    }
    void prepare(SoGroup &parent, bool append, bool prepareChildren = true)
    {
	this->hierarchy.prepare(prepareChildren,
	    this->effect == Effect::Structural ? PeerEffects::Revision::Structural :
		PeerEffects::Revision::Frame);
	auto &state = *this->scene.d;
	bool needsSceneWalk = !append;
	auto collect = [&](SoNode *node, SoGroup *) {
	    if (node->isOfType(SoBRLSceneGroup::getClassTypeId())) {
		const auto key = normalized_index_key(static_cast<SoBRLSceneGroup *>(node)->groupPath.getValue().getString());
		this->groups.set(key, nullptr);
		needsSceneWalk = needsSceneWalk || state.indexes.groupPathIndex.count(key);
	    }
	    if (!node->isOfType(SoBRLDatabaseSource::getClassTypeId())) return;
	    auto *source = static_cast<SoBRLDatabaseSource *>(node);
	    if (!this->affected.insert(source).second) return;
	    this->sources.push_back(source);
	    const auto path = normalized_index_key(source->path.getValue().getString());
	    const auto key = normalized_index_key(database_source_effective_instance_key(source).getString());
	    this->paths.set(path, nullptr); this->instances.set(key, nullptr); this->parents.set(key, nullptr);
	    this->pathValues.emplace(path, std::vector<SoBRLDatabaseSource *>());
	    needsSceneWalk = needsSceneWalk || state.indexes.databaseSourcePathInstancesIndex.count(path) ||
		state.indexes.databaseSourceInstanceIndex.count(key) || state.indexes.databaseSourceOrderIndex.count(source);
	};
	for (const auto &child : this->retainedChildren)
	    if (!this->walk(child.get(), &parent, collect)) return;

	// Only colliding keys and removals need a reachability walk. Ordinary
	// insertion prepares storage proportional to the incoming subtree.
	if (this->sources.empty() && this->groups.empty()) needsSceneWalk = false;
	std::unordered_set<SoBRLDatabaseSource *> listed;
	auto resolve = [&](SoNode *node, SoGroup *nodeParent) {
	    if (node->isOfType(SoBRLSceneGroup::getClassTypeId())) {
		const auto key = normalized_index_key(static_cast<SoBRLSceneGroup *>(node)->groupPath.getValue().getString());
		if (this->groups.contains(key)) this->groups.set(key, static_cast<SoGroup *>(node));
	    }
	    if (!node->isOfType(SoBRLDatabaseSource::getClassTypeId())) return;
	    auto *source = static_cast<SoBRLDatabaseSource *>(node);
	    if (this->affected.count(source)) this->reachable.insert(source);
	    const auto path = normalized_index_key(source->path.getValue().getString());
	    const auto key = normalized_index_key(database_source_effective_instance_key(source).getString());
	    auto values = this->pathValues.find(path);
	    if (values != this->pathValues.end()) {
		this->paths.set(path, source);
		if (listed.insert(source).second) values->second.push_back(source);
	    }
	    if (this->instances.contains(key)) {
		this->instances.set(key, source); this->parents.set(key, nodeParent);
	    }
	};
	if (needsSceneWalk) {
	    if (!this->walk(state.root, nullptr, resolve)) return;
	} else {
	    for (const auto &child : this->retainedChildren)
		if (!this->walk(child.get(), &parent, resolve)) return;
	}
	for (auto *source : this->sources) {
	    if (!this->reachable.count(source)) this->retired.push_back(source);
	    else if (!state.indexes.databaseSourceOrderIndex.count(source)) {
		this->added.push_back(source);
		this->newOrder.emplace(source, 0);
		this->newRouting.emplace(source->getCompactSourceRoutingId(), source);
	    }
	}
	if (!this->added.empty()) {
	    state.indexes.databaseSourceOrder.reserve(state.indexes.databaseSourceOrder.size() + this->added.size());
	    state.indexes.databaseSourceOrderIndex.reserve(state.indexes.databaseSourceOrderIndex.size() + this->added.size());
	    state.indexes.databaseSourceRoutingIndex.reserve(state.indexes.databaseSourceRoutingIndex.size() + this->added.size());
	}
	size_t newPaths = 0;
	for (const auto &entry : this->pathValues)
	    if (!entry.second.empty() && !state.indexes.databaseSourcePathInstancesIndex.count(entry.first)) ++newPaths;
	if (newPaths) state.indexes.databaseSourcePathInstancesIndex.reserve(state.indexes.databaseSourcePathInstancesIndex.size() + newPaths);
	this->groups.prepare(); this->paths.prepare(); this->instances.prepare(); this->parents.prepare();
	this->valid = true;
    }
    template <typename Visit>
    bool walk(SoNode *root, SoGroup *parent, Visit visit) const
    {
	return scene_walk(root, parent, visit, [this](SoGroup *group) { return this->hierarchy.preparedChildren(group); });
    }
    BObolSceneController &scene;
    std::vector<SbModernUtils::SoNodeRef> retainedChildren;
    HierarchyPublication hierarchy;
    PreparedSceneIndex<SoGroup> groups;
    PreparedSceneIndex<SoBRLDatabaseSource> paths, instances;
    PreparedSceneIndex<SoGroup> parents;
    decltype(SceneIndexes::databaseSourcePathInstancesIndex) pathValues;
    decltype(SceneIndexes::databaseSourceOrderIndex) newOrder;
    decltype(SceneIndexes::databaseSourceRoutingIndex) newRouting;
    std::unordered_set<SoBRLDatabaseSource *> affected, reachable;
    std::vector<SoBRLDatabaseSource *> sources, added, retired;
    bool valid = false;
    Effect effect = Effect::Structural;
};

class BObolSceneController::SourceChildEffects : public BObolSourceChildEffects {
public:
    explicit SourceChildEffects(BObolSceneController &owner) : scene(owner) {}
    void stageChildOrder(SoGroup &parent, const std::vector<SoNode *> &children) override
    {
	this->orders[&parent] = children;
    }
    void stageFrameEffect(SoNode &node) override
    {
	if (std::none_of(this->frameNodes.begin(), this->frameNodes.end(),
	    [&node](const auto &current) { return current.get() == &node; }))
	    this->frameNodes.emplace_back(&node);
    }
    void prepare() override
    {
	this->publication = std::make_unique<ChildPublication>(this->scene, this->orders);
	for (const auto &entry : this->orders)
	    this->publication->includeStructuralEffect(entry.first);
	for (const auto &node : this->frameNodes) this->publication->includeEffect(node.get());
    }
    void commit() noexcept override { this->publication->commitEffects(); }
private:
    BObolSceneController &scene;
    SceneChildOrders orders;
    std::vector<SbModernUtils::SoNodeRef> frameNodes;
    std::unique_ptr<ChildPublication> publication;
};

std::unique_ptr<BObolSourceRealizationEffects>
BObolSceneController::prepareRealizationEffects(
    BObolSceneController &scene, BObolSourceRealizationEffects *downstream)
{
    class Effects : public BObolSourceRealizationEffects {
    public:
	Effects(BObolSceneController &owner, BObolSourceRealizationEffects *next) :
	    scene(owner), downstream(next) {}
	void prepare(const SoBRLDatabaseSource &source, bool realized, const SbString &diagnostic,
	    const std::vector<SoNode *> &children, const std::vector<SoNode *> &changedNodes) override
	{
	    if (realized) {
		this->peers.reset();
		this->childEffects = std::make_unique<ChildPublication>(this->scene,
		    const_cast<SoBRLDatabaseSource &>(source), children,
		    ChildPublication::Effect::Frame);
		for (auto *node : changedNodes) this->childEffects->includeEffect(node);
	    } else {
		this->childEffects.reset();
		this->peers = std::make_unique<PeerEffects>(
		    this->scene, const_cast<SoBRLDatabaseSource *>(&source));
	    }
	    if (this->downstream)
		this->downstream->prepare(source, realized, diagnostic, children, changedNodes);
	}
	void commit(bool changed) noexcept override
	{
	    if (changed) {
		if (this->childEffects) this->childEffects->commitEffects();
		else PeerEffects::frameCommitted(this->peers.get());
	    }
	    if (this->downstream) this->downstream->commit(changed);
	}
	void notify() override { if (this->downstream) this->downstream->notify(); }
    private:
	BObolSceneController &scene;
	BObolSourceRealizationEffects *downstream;
	std::unique_ptr<PeerEffects> peers;
	std::unique_ptr<ChildPublication> childEffects;
    };
    return std::make_unique<Effects>(scene, downstream);
}

int
BObolSceneController::appendChildToGroup(const char *groupPath,
	SoNode *child)
{
    if (!child)
	return -1;

    SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group)
	return -1;
    if (scene_group_find_child_index(group, child) >= 0)
	return 0;

    ChildPublication publication(*this, *group, *child, true);
    return publication.publish();
}

int
BObolSceneController::removeChildFromGroup(const char *groupPath,
	SoNode *child)
{
    if (!child)
	return -1;

    SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group)
	return -1;

    const int childIndex = scene_group_find_child_index(group, child);
    if (childIndex < 0)
	return 0;

    ChildPublication publication(*this, *group, *child, false);
    return publication.publish();
}

int
BObolSceneController::eraseGroupSubpath(const char *parentGroupPath,
					const char *subpath)
{
    SoGroup *parent = this->findIndexedGroup(parentGroupPath);
    if (!parent)
	return -1;

    return this->removeGroupSubpath(parent, subpath);
}

int
BObolSceneController::removeGroup(const char *groupPath)
{
    return this->removeGroupSubpath(scene_root_group(this->d->root), groupPath);
}

int
BObolSceneController::removeGroupSubpath(SoGroup *parent, const char *subpath)
{
    std::vector<std::string> components;
    scene_group_path_components(subpath, components);
    if (components.empty())
	return 0;

    if (!parent)
	return -1;

    for (size_t i = 0; i + 1 < components.size(); i++) {
	parent = scene_group_find_child(parent, components[i].c_str());
	if (!parent)
	    return 0;
    }

    SoGroup *target = scene_group_find_child(parent,
		      components.back().c_str());
    if (!target)
	return 0;

    ChildPublication publication(*this, *parent, *target, false);
    return publication.publish();
}

int
BObolSceneController::clearGroup(const char *groupPath)
{
    SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group)
	return -1;

    const int removed = group->getNumChildren();
    if (removed <= 0)
	return 0;

    ChildPublication publication(*this, *group);
    return publication.publish() < 0 ? -1 : removed;
}

int
BObolSceneController::getGroupChildCount(const char *groupPath) const
{
    const SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group)
	return -1;
    return group->getNumChildren();
}

int
BObolSceneController::getGroupDescendantGroupCount(
    const char *groupPath) const
{
    const SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group)
	return -1;
    return count_scene_groups_recursive(group);
}

int
BObolSceneController::getGroupDatabaseSourceCount(
    const char *groupPath) const
{
    const SoGroup *group = this->findIndexedGroup(groupPath);
    if (!group)
	return -1;
    return count_database_sources_recursive(group);
}

SoNode *
BObolSceneController::findShape(const char *shapePath) const
{
    return scene_shape_find_path(this->d->root, shapePath, NULL);
}

SoGroup *
BObolSceneController::findShapeParent(const char *shapePath) const
{
    SoGroup *parent = NULL;
    (void)scene_shape_find_path(this->d->root, shapePath, &parent);
    return parent;
}

int
BObolSceneController::moveShapeToGroup(const char *shapePath,
				       const char *groupPath)
{
    SoGroup *currentParent = NULL;
    SoNode *shape = scene_shape_find_path(this->d->root, shapePath,
					  &currentParent);
    if (!shape)
	return 0;

    SoGroup *targetGroup = this->findIndexedGroup(groupPath);
    if (!targetGroup)
	return -1;
    if (targetGroup == currentParent)
	return 0;
    if (scene_group_find_child_index(targetGroup, shape) >= 0)
	return 0;

    const int currentIndex =
	scene_group_find_child_index(currentParent, shape);
    if (currentIndex < 0)
	return -1;

    return ChildPublication::moveShape(*this, *currentParent, *targetGroup, *shape);
}

int
BObolSceneController::removeShape(const char *shapePath)
{
    SoGroup *parent = NULL;
    SoNode *shape = scene_shape_find_path(this->d->root, shapePath, &parent);
    if (!shape)
	return 0;

    const int childIndex = scene_group_find_child_index(parent, shape);
    if (childIndex < 0)
	return -1;

    ChildPublication publication(*this, *parent, *shape, false);
    return publication.publish();
}

class BObolSceneController::ShapePublication {
public:
    ShapePublication(BObolSceneController &owner, SoNode &target) : scene(owner), retained(&target) {}

    template <typename Field, typename Value>
    void set(Field &field, const Value &requested)
    {
	if constexpr (std::is_same_v<Field, SoSFFloat>) {
	    if (!scene_group_float_different(field.getValue(), requested)) return;
	} else if (field.getValue() == requested) return;
	auto next = std::make_unique<Field>();
	next->enableNotify(FALSE);
	next->setValue(requested);
	const SoField *value = next.get();
	this->owned.push_back(std::move(next));
	if (this->values.prepare(field, value)) this->changes.push_back({&field, true});
    }
    int publish()
    {
	if (this->changes.empty()) return 0;
	PeerEffects effects(this->scene, this->retained.get());
	PreparedNotifications<std::vector<PublicationFieldChange>> notifications(*this->retained.get(), std::move(this->changes));
	this->values.commit();
	effects.commit();
	this->scene.advanceFrameRevision();
	notifications.restore();
	std::exception_ptr failure;
	notifications.notify(failure);
	if (failure) std::rethrow_exception(failure);
	return 1;
    }
private:
    BObolSceneController &scene;
    SbModernUtils::SoNodeRef retained;
    // Stage only requested fields: geometry and unrelated derived metadata
    // belong to their existing owners and must retain their identity.
    std::vector<std::unique_ptr<SoField>> owned;
    PreparedScalarValues values;
    std::vector<PublicationFieldChange> changes;
};

template <typename Configure>
int
BObolSceneController::publishShapeState(const char *shapePath, Configure configure)
{
    SoNode *shape = scene_shape_find_path(this->d->root, shapePath, nullptr);
    if (!shape) return -1;
    ShapePublication publication(*this, *shape);
    if (shape->isOfType(SoBRLVListShape::getClassTypeId()))
	configure(*static_cast<SoBRLVListShape *>(shape), publication);
    else
	configure(*static_cast<SoBRLMeshShape *>(shape), publication);
    return publication.publish();
}

int
BObolSceneController::setShapeDrawState(const char *shapePath,
					int drawMode,
					SbBool databaseIntent,
					SbBool overlayIntent,
					SbBool hudIntent)
{
    return this->publishShapeState(shapePath, [&](auto &shape, auto &publication) {
	publication.set(shape.drawMode, drawMode);
	publication.set(shape.databaseIntent, databaseIntent);
	publication.set(shape.overlayIntent, overlayIntent);
	publication.set(shape.hudIntent, hudIntent);
    });
}

int
BObolSceneController::setShapeDisplayState(const char *shapePath,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision)
{
    return this->publishShapeState(shapePath, [&](auto &shape, auto &publication) {
	publication.set(shape.visible, visible);
	publication.set(shape.selected, selected);
	publication.set(shape.highlighted, highlighted);
	publication.set(shape.lineStyle, lineStyle);
	publication.set(shape.lineWidth, lineWidth);
	publication.set(shape.transparency, transparency);
	publication.set(shape.colorOverride, colorOverride);
	if (!scene_group_color_equal(shape.color.getValue(), color)) publication.set(shape.color, color);
	publication.set(shape.materialColorValid, materialColorValid);
	if (!scene_group_color_equal(shape.materialColor.getValue(), materialColor)) publication.set(shape.materialColor, materialColor);
	publication.set(shape.materialRevision, materialRevision);
    });
}

int
BObolSceneController::setShapePlacementState(const char *shapePath,
	SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize)
{
    return this->publishShapeState(shapePath, [&](auto &shape, auto &publication) {
	publication.set(shape.drawMatrixValid, drawMatrixValid);
	if (!shape.drawMatrix.getValue().equals(drawMatrix, scene_state_tolerance)) publication.set(shape.drawMatrix, drawMatrix);
	publication.set(shape.drawCenterValid, drawCenterValid);
	if (!scene_group_vec3f_equal(shape.drawCenter.getValue(), drawCenter)) publication.set(shape.drawCenter, drawCenter);
	publication.set(shape.drawSizeValid, drawSizeValid);
	publication.set(shape.drawSize, drawSize);
    });
}

int
BObolSceneController::publishDatabaseSourceAuxiliaryLineSet(
    const char *sourcePath,
    const char *name,
    const SbVec3f *points,
    const int32_t *commands,
    int count,
    const BObolAuxiliaryLineSetDisplayState *displayState)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->publishDatabaseSourceInstanceAuxiliaryLineSet(
	       sourceInstanceKey.getString(),
	       name, points, commands, count, displayState);
}

int
BObolSceneController::publishDatabaseSourceAuxiliarySourceLineSet(
    const char *sourcePath,
    const char *auxiliarySourcePath,
    const char *displayName,
    const SbVec3f *points,
    const int32_t *commands,
    int count,
    const BObolAuxiliaryLineSetDisplayState *displayState)
{
    return this->publishDatabaseSourceInstanceAuxiliarySourceLineSet(
	       sourcePath, auxiliarySourcePath, displayName, points, commands,
	       count, displayState);
}

int
BObolSceneController::publishDatabaseSourceInstanceAuxiliaryLineSet(
    const char *sourceInstanceKey,
    const char *name,
    const SbVec3f *points,
    const int32_t *commands,
    int count,
    const BObolAuxiliaryLineSetDisplayState *displayState)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0] || !name || !name[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    SourceChildEffects effects(*this);
    return source->publishAuxiliaryLineSet(name, points, commands, count, displayState, &effects);
}

int
BObolSceneController::publishDatabaseSourceInstanceAuxiliarySourceLineSet(
    const char *sourceInstanceKey,
    const char *auxiliarySourcePath,
    const char *displayName,
    const SbVec3f *points,
    const int32_t *commands,
    int count,
    const BObolAuxiliaryLineSetDisplayState *displayState)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0] ||
	!auxiliarySourcePath || !auxiliarySourcePath[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    SourceChildEffects effects(*this);
    return source->publishAuxiliarySourceLineSet(auxiliarySourcePath, displayName,
	points, commands, count, displayState, &effects);
}

int
BObolSceneController::publishDatabaseSourceExternalLineSet(
    const char *sourcePath,
    const BObolExternalLineSet &lineSet)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->publishDatabaseSourceInstanceExternalLineSet(
	       sourceInstanceKey.getString(), lineSet);
}

int
BObolSceneController::publishDatabaseSourceInstanceExternalLineSet(
    const char *sourceInstanceKey,
    const BObolExternalLineSet &lineSet)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    SourceChildEffects effects(*this);
    return source->publishExternalLineSet(lineSet, &effects);
}

int
BObolSceneController::publishDatabaseSourceInstancePrimitiveWireframe(
    const char *sourceInstanceKey,
    struct rt_db_internal *intern,
    const struct bg_tess_tol *ttol,
    const struct bn_tol *tol)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0] || !intern)
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    SourceChildEffects effects(*this);
    return source->publishPrimitiveWireframe(intern, ttol, tol, &effects);
}

int
BObolSceneController::publishDatabaseSourceExternalPointSet(
    const char *sourcePath,
    const BObolExternalPointSet &pointSet)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->publishDatabaseSourceInstanceExternalPointSet(
	       sourceInstanceKey.getString(), pointSet);
}

int
BObolSceneController::publishDatabaseSourceInstanceExternalPointSet(
    const char *sourceInstanceKey,
    const BObolExternalPointSet &pointSet)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    SourceChildEffects effects(*this);
    return source->publishExternalPointSet(pointSet, &effects);
}

int
BObolSceneController::publishDatabaseSourceExternalTriangleMesh(
    const char *sourcePath,
    const BObolExternalTriangleMesh &triangleMesh)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->publishDatabaseSourceInstanceExternalTriangleMesh(
	       sourceInstanceKey.getString(), triangleMesh);
}

int
BObolSceneController::publishDatabaseSourceInstanceExternalTriangleMesh(
    const char *sourceInstanceKey,
    const BObolExternalTriangleMesh &triangleMesh)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    SourceChildEffects effects(*this);
    return source->publishExternalTriangleMesh(triangleMesh, &effects);
}

int
BObolSceneController::publishDatabaseSourceExternalAnnotation(
    const char *sourcePath,
    const BObolExternalAnnotation &annotation)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->publishDatabaseSourceInstanceExternalAnnotation(
	       sourceInstanceKey.getString(), annotation);
}

int
BObolSceneController::publishDatabaseSourceInstanceExternalAnnotation(
    const char *sourceInstanceKey,
    const BObolExternalAnnotation &annotation)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    SourceChildEffects effects(*this);
    return source->publishExternalAnnotation(annotation, &effects);
}

int
BObolSceneController::clearDatabaseSourceExternalPrimaryGeometry(
    const char *sourcePath)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->clearDatabaseSourceInstanceExternalPrimaryGeometry(
	       sourceInstanceKey.getString());
}

int
BObolSceneController::clearDatabaseSourceInstanceExternalPrimaryGeometry(
    const char *sourceInstanceKey)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    SourceChildEffects effects(*this);
    return source->clearExternalPrimaryGeometry(&effects);
}

int
BObolSceneController::clearDatabaseSourceAuxiliaryShapes(const char *sourcePath)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->clearDatabaseSourceInstanceAuxiliaryShapes(
	       sourceInstanceKey.getString());
}

int
BObolSceneController::clearDatabaseSourceInstanceAuxiliaryShapes(
    const char *sourceInstanceKey)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return -1;

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    SourceChildEffects effects(*this);
    return source->removeAuxiliaryShapes(&effects);
}

int
BObolSceneController::setShapeSourceState(const char *shapePath,
	const char *ownerSourcePath,
	uint32_t ownerSourceRevision,
	uint32_t ownerInputsRevision,
	uint32_t ownerViewRevision,
	uint32_t ownerRealizedRevision,
	uint32_t ownerRealizedSourceRevision,
	uint32_t ownerRealizedInputsRevision,
	uint32_t ownerRealizedViewRevision,
	int ownerRealizationStatus,
	const char *ownerRealizationDiagnostic,
	const char *ownerRealizationIdentity,
	SbBool ownerSourceStale,
	uint32_t ownerStaleReason)
{
    return this->publishShapeState(shapePath, [&](auto &shape, auto &publication) {
	publication.set(shape.ownerSourcePath, ownerSourcePath ? ownerSourcePath : "");
	publication.set(shape.ownerSourceRevision, ownerSourceRevision);
	publication.set(shape.ownerInputsRevision, ownerInputsRevision);
	publication.set(shape.ownerViewRevision, ownerViewRevision);
	publication.set(shape.ownerRealizedRevision, ownerRealizedRevision);
	publication.set(shape.ownerRealizedSourceRevision, ownerRealizedSourceRevision);
	publication.set(shape.ownerRealizedInputsRevision, ownerRealizedInputsRevision);
	publication.set(shape.ownerRealizedViewRevision, ownerRealizedViewRevision);
	publication.set(shape.ownerRealizationStatus, ownerRealizationStatus);
	publication.set(shape.ownerRealizationDiagnostic, ownerRealizationDiagnostic ? ownerRealizationDiagnostic : "");
	publication.set(shape.ownerRealizationIdentity, ownerRealizationIdentity ? ownerRealizationIdentity : "");
	publication.set(shape.ownerSourceStale, ownerSourceStale);
	publication.set(shape.ownerStaleReason, ownerStaleReason);
    });
}

SoBRLDatabaseSource *
BObolSceneController::getDatabaseSource(int index) const
{
    if (index < 0 || !this->d->root ||
	!this->d->root->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    if (!this->d->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();
    const size_t pos = static_cast<size_t>(index);
    if (pos >= this->d->indexes.databaseSourceOrder.size())
	return NULL;
    return this->d->indexes.databaseSourceOrder[pos];
}

int
BObolSceneController::getDatabaseSourceCount(void) const
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return 0;

    if (!this->d->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();
    return static_cast<int>(this->d->indexes.databaseSourceOrder.size());
}

SoBRLDatabaseSource *
BObolSceneController::findDatabaseSourceRoutingId(uint64_t routingId) const
{
    if (!routingId || !this->d->root ||
	!this->d->root->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    if (!this->d->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();
    const auto route = this->d->indexes.databaseSourceRoutingIndex.find(routingId);
    return route != this->d->indexes.databaseSourceRoutingIndex.end() ?
	route->second : NULL;
}

SoBRLDatabaseSource *
BObolSceneController::findDatabaseSource(const char *sourcePath) const
{
    if (!sourcePath || !sourcePath[0] || !this->d->root ||
	!this->d->root->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    return this->findIndexedDatabaseSource(sourcePath);
}

SoBRLDatabaseSource *
BObolSceneController::findDatabaseSourceInstance(
    const char *sourceInstanceKey) const
{
    if (!sourceInstanceKey || !sourceInstanceKey[0] || !this->d->root ||
	!this->d->root->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    return this->findIndexedDatabaseSourceInstance(sourceInstanceKey);
}

int
BObolSceneController::replaceDatabaseSource(const char *sourcePath,
	struct db_i *database,
	int drawMode,
	uint32_t sourceRevision)
{
    return this->replaceDatabaseSourceInstance(sourcePath, sourcePath,
	    database, drawMode, sourceRevision);
}

int
BObolSceneController::replaceDatabaseSourceInstance(
    const char *sourceInstanceKey,
    const char *sourcePath,
    struct db_i *database,
    int drawMode,
    uint32_t sourceRevision)
{
    return this->replaceDatabaseSourceInstanceRepresentation(
	       sourceInstanceKey, sourcePath, NULL, -1, database, drawMode,
	       sourceRevision);
}

int
BObolSceneController::replaceDatabaseSourceInstanceRepresentation(
    const char *sourceInstanceKey,
    const char *sourcePath,
    const char *sourceRepresentationKey,
    int sourceRepresentationMode,
    struct db_i *database,
    int drawMode,
    uint32_t sourceRevision)
{
    BObolDatabaseSourcePublishState state;
    state.sourceInstanceKey = sourceInstanceKey;
    state.sourcePath = sourcePath;
    state.sourceRepresentationKey = sourceRepresentationKey;
    state.database = database;
    state.drawMode = drawMode;
    state.representationMode = sourceRepresentationMode;
    state.sourceRevisionValid = sourceRevision ? TRUE : FALSE;
    state.sourceRevision = sourceRevision;
    return this->publishDatabaseSourceInstance(state);
}

int
BObolSceneController::renameDatabaseSource(const char *sourcePath,
	const char *newSourcePath,
	uint32_t sourceRevision)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoBRLDatabaseSource *source = this->findIndexedDatabaseSource(sourcePath);
    SbString sourceInstanceKey = database_source_effective_instance_key(source);
    if (sourceInstanceKey.getLength() == 0)
	return 0;

    return this->renameDatabaseSourceInstance(sourceInstanceKey.getString(),
	    newSourcePath,
	    newSourcePath, sourceRevision);
}

class BObolSceneController::SourceIndexPublication {
public:
    SourceIndexPublication(BObolSceneController &owner, HierarchyPublication &preparedHierarchy, SoBRLDatabaseSource &renamed, SoGroup &parent,
	const char *nextInstance, const char *nextPath, SoBRLDatabaseSource *conflict,
	SoGroup *conflictParent, int conflictIndex) :
	scene(owner), source(renamed), retainedSource(&renamed), paths(owner.d->indexes.databaseSourcePathIndex),
	instances(owner.d->indexes.databaseSourceInstanceIndex), parents(owner.d->indexes.databaseSourceInstanceParentIndex),
	groups(owner.d->indexes.groupPathIndex), hierarchy(preparedHierarchy)
    {
	auto &state = *owner.d;
	if (conflict) {
	    std::unordered_set<SoNode *> unreachable;
	    std::unordered_map<SoBRLDatabaseSource *, SoGroup *> survivingParents;
	    collectSubtree(conflict, unreachable);
	    removeReachable(state.root, nullptr, conflictParent, conflictIndex, unreachable, survivingParents);
	    // Replacing this edge must retire the conflicting instance without
	    // removing the source being renamed or leaving a duplicate identity.
	    if (unreachable.count(&renamed) || !unreachable.count(conflict)) return;
	    for (const auto &entry : survivingParents)
		this->parents.set(normalized_index_key(database_source_effective_instance_key(entry.first).getString()), entry.second);
	    for (auto *candidate : state.indexes.databaseSourceOrder)
		if (unreachable.count(candidate)) this->retired.push_back(candidate);
	    for (const auto &entry : state.indexes.groupPathIndex)
		if (unreachable.count(entry.second)) this->groups.set(entry.first, nullptr);
	    this->hierarchy.removeChild(conflictParent, conflictIndex);
	}
	this->unindexed.insert(&renamed);
	this->unindexed.insert(this->retired.begin(), this->retired.end());
	for (auto *candidate : this->unindexed) {
	    this->pathChanges.emplace(normalized_index_key(candidate->path.getValue().getString()), false);
	    const auto key = normalized_index_key(database_source_effective_instance_key(candidate).getString());
	    this->instances.set(key, nullptr);
	    this->parents.set(key, nullptr);
	}
	const auto key = normalized_index_key(nextInstance);
	this->instances.set(key, &renamed);
	this->parents.set(key, &parent);
	this->pathChanges[normalized_index_key(nextPath)] = true;
	for (const auto &entry : this->pathChanges) {
	    auto found = state.indexes.databaseSourcePathInstancesIndex.find(entry.first);
	    SoBRLDatabaseSource *last = nullptr;
	    if (found != state.indexes.databaseSourcePathInstancesIndex.end()) {
		size_t nextCount = entry.second ? 1 : 0;
		for (auto *candidate : found->second)
		    if (!this->unindexed.count(candidate)) { last = candidate; ++nextCount; }
		if (entry.second) found->second.reserve(nextCount);
	    } else if (entry.second) this->newPaths[entry.first].push_back(&renamed);
	    this->paths.set(entry.first, entry.second ? &renamed : last);
	}
	if (!this->newPaths.empty())
	    state.indexes.databaseSourcePathInstancesIndex.reserve(state.indexes.databaseSourcePathInstancesIndex.size() + this->newPaths.size());
	if (state.indexes.databaseSourceOrderIndex.find(&renamed) == state.indexes.databaseSourceOrderIndex.end()) {
	    this->newOrder.emplace(&renamed, state.indexes.databaseSourceOrder.size());
	    state.indexes.databaseSourceOrderIndex.reserve(state.indexes.databaseSourceOrderIndex.size() + 1);
	    state.indexes.databaseSourceOrder.reserve(state.indexes.databaseSourceOrder.size() + 1);
	}
	const auto route = renamed.getCompactSourceRoutingId();
	if (state.indexes.databaseSourceRoutingIndex.find(route) == state.indexes.databaseSourceRoutingIndex.end()) {
	    this->newRouting.emplace(route, &renamed);
	    state.indexes.databaseSourceRoutingIndex.reserve(state.indexes.databaseSourceRoutingIndex.size() + 1);
	}
	this->paths.prepare(); this->instances.prepare(); this->parents.prepare(); this->groups.prepare();
	this->hierarchy.effects().include(&renamed, PeerEffects::Revision::Structural);
	this->valid = true;
    }
    bool isValid() const { return this->valid; }
    static void committed(void *context) noexcept
    {
	auto &self = *static_cast<SourceIndexPublication *>(context);
	auto &state = *self.scene.d;
	self.hierarchy.commit();
	for (const auto &entry : self.pathChanges) {
	    auto found = state.indexes.databaseSourcePathInstancesIndex.find(entry.first);
	    if (found == state.indexes.databaseSourcePathInstancesIndex.end()) {
		if (entry.second) state.indexes.databaseSourcePathInstancesIndex.insert(self.newPaths.extract(entry.first));
		continue;
	    }
	    auto &values = found->second;
	    values.erase(std::remove_if(values.begin(), values.end(), [&self](SoBRLDatabaseSource *candidate) {
		return self.unindexed.count(candidate) != 0;
	    }), values.end());
	    if (entry.second) values.push_back(&self.source);
	    if (values.empty()) state.indexes.databaseSourcePathInstancesIndex.erase(found);
	}
	self.paths.commit(); self.instances.commit(); self.parents.commit(); self.groups.commit();
	for (auto *retired : self.retired) {
	    state.indexes.databaseSourceRoutingIndex.erase(retired->getCompactSourceRoutingId());
	    state.indexes.databaseSourceOrderIndex.erase(retired);
	}
	if (!self.newOrder.empty()) state.indexes.databaseSourceOrderIndex.insert(self.newOrder.extract(self.newOrder.begin()));
	if (!self.newRouting.empty()) state.indexes.databaseSourceRoutingIndex.insert(self.newRouting.extract(self.newRouting.begin()));
	auto &order = state.indexes.databaseSourceOrder;
	order.erase(std::remove_if(order.begin(), order.end(), [&self](SoBRLDatabaseSource *candidate) {
	    return self.unindexed.count(candidate) != 0;
	}), order.end());
	order.push_back(&self.source);
	for (size_t i = 0; i < order.size(); ++i) state.indexes.databaseSourceOrderIndex.find(order[i])->second = i;
	self.scene.advanceStructuralRevision();
    }
    void notify(std::exception_ptr &failure)
    {
	this->hierarchy.notify(failure);
    }
    static void collectSubtree(SoNode *node, std::unordered_set<SoNode *> &nodes)
    {
	if (!node || !nodes.insert(node).second) return;
	if (node->isOfType(SoGroup::getClassTypeId())) {
	    auto *group = static_cast<SoGroup *>(node);
	    for (int i = 0; i < group->getNumChildren(); ++i) collectSubtree(group->getChild(i), nodes);
	}
    }
private:
    static void removeReachable(SoNode *node, SoGroup *parent, SoGroup *removedParent, int removedIndex,
	std::unordered_set<SoNode *> &unreachable, std::unordered_map<SoBRLDatabaseSource *, SoGroup *> &survivingParents)
    {
	if (!node) return;
	if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    auto *source = static_cast<SoBRLDatabaseSource *>(node);
	    if (unreachable.count(node) || survivingParents.count(source)) survivingParents[source] = parent;
	}
	unreachable.erase(node);
	if (node->isOfType(SoGroup::getClassTypeId())) {
	    auto *group = static_cast<SoGroup *>(node);
	    for (int i = 0; i < group->getNumChildren(); ++i)
		if (group != removedParent || i != removedIndex)
		    removeReachable(group->getChild(i), group, removedParent, removedIndex, unreachable, survivingParents);
	}
    }
    BObolSceneController &scene;
    SoBRLDatabaseSource &source;
    SbModernUtils::SoNodeRef retainedSource;
    PreparedSceneIndex<SoBRLDatabaseSource> paths, instances;
    PreparedSceneIndex<SoGroup> parents, groups;
    std::vector<SoBRLDatabaseSource *> retired;
    std::unordered_set<SoBRLDatabaseSource *> unindexed;
    std::unordered_map<std::string, bool> pathChanges;
    decltype(SceneIndexes::databaseSourcePathInstancesIndex) newPaths;
    decltype(SceneIndexes::databaseSourceOrderIndex) newOrder;
    decltype(SceneIndexes::databaseSourceRoutingIndex) newRouting;
    HierarchyPublication &hierarchy;
    bool valid = false;
};

class BObolSceneController::SourcePublication {
public:
    SourcePublication(BObolSceneController &owner, const char *currentSourceInstanceKey,
	bool requireExisting, const BObolDatabaseSourcePublishState &state,
	const BObolSceneGroupPublishState *groupState) :
	scene(owner), key(state.sourceInstanceKey), path(state.sourcePath), sourceOwner(nullptr),
	hierarchy(owner)
    {
	const char *currentKey = currentSourceInstanceKey &&
	    currentSourceInstanceKey[0] ? currentSourceInstanceKey :
	    this->key.c_str();
	this->source = owner.findIndexedDatabaseSourceInstance(currentKey);
	if (requireExisting && !this->source) return;
	SoBRLDatabaseSource *destination =
	    owner.findIndexedDatabaseSourceInstance(this->key.c_str());
	if (destination && destination != this->source) return;
	SoGroup *parent = this->source ?
	    owner.findIndexedDatabaseSourceInstanceParent(currentKey) : nullptr;
	const bool insertion = !this->source;
	if (insertion) this->source = new SoBRLDatabaseSource;
	this->sourceOwner = SbModernUtils::SoNodeRef(this->source);
	if (!insertion && !parent) return;
	this->prepare(parent, state.targetGroupPath, groupState);
    }
    SourcePublication(BObolSceneController &owner, SoBRLDatabaseSource &candidate,
	SoGroup &parent, const char *targetGroupPath) :
	scene(owner), key(database_source_effective_instance_key(&candidate).getString()),
	path(candidate.path.getValue().getString()), sourceOwner(&candidate), source(&candidate), hierarchy(owner)
    {
	this->prepare(&parent, targetGroupPath[0] ? targetGroupPath : "/", nullptr);
    }
    bool isValid() const { return this->valid; }
    SoBRLDatabaseSource &getSource() const { return *this->source; }
    bool changesScene() const { return this->structural || (this->metadata && this->metadata->changed()); }
    bool didCommit() const { return this->published; }
    static void committed(void *context) noexcept
    {
	auto &self = *static_cast<SourcePublication *>(context);
	if (self.metadata) self.metadata->commit();
	if (self.indexes) SourceIndexPublication::committed(self.indexes.get());
	else {
	    self.hierarchy.commit();
	    self.scene.advanceFrameRevision();
	}
	self.published = true;
	// Source observers may immediately edit this group. Restore its fields
	// before returning to the source's notification phase.
	if (self.metadata) self.metadata->restore();
    }
    void notify(std::exception_ptr &failure)
    {
	if (!this->published) return;
	if (this->metadata) this->metadata->notify(failure);
	this->hierarchy.notify(failure);
    }
private:
    void prepare(SoGroup *parent, const char *targetGroupPath, const BObolSceneGroupPublishState *groupState)
    {
	if (parent && scene_group_find_child_index(parent, this->source) < 0) return;
	const bool insertion = !parent;
	auto *root = static_cast<SoGroup *>(this->scene.d->root);
	SoGroup *target = parent ? parent : root;
	if (targetGroupPath && targetGroupPath[0])
	    target = this->hierarchy.ensurePath(targetGroupPath);
	if (parent && parent != target) {
	    if (this->hierarchy.reaches(this->source, target)) return;
	    auto &previous = this->hierarchy.nextChildren(parent);
	    previous.erase(std::find(previous.begin(), previous.end(), this->source));
	}
	if (insertion || parent != target) this->hierarchy.nextChildren(target).push_back(this->source);
	this->structural = insertion || parent != target ||
	    !database_source_instance_key_equal(this->source, this->key.c_str()) ||
	    !database_source_path_equal(this->source, this->path.c_str());
	if (this->structural) {
	    this->indexes = std::make_unique<SourceIndexPublication>(this->scene, this->hierarchy, *this->source, *target,
		this->key.c_str(), this->path.c_str(), nullptr, nullptr, -1);
	    if (!this->indexes->isValid()) return;
	} else this->hierarchy.effects().include(this->source, PeerEffects::Revision::Frame);
	if (groupState) {
	    if (!target->isOfType(SoBRLSceneGroup::getClassTypeId())) return;
	    this->metadata = std::make_unique<GroupPublication>(this->scene, *static_cast<SoBRLSceneGroup *>(target));
	    this->metadata->setIntent(targetGroupPath, groupState->intentPath, groupState->drawMode,
		groupState->fallbackDrawMode, groupState->overlayIntent, groupState->revalidationRevision);
	    this->metadata->setDisplay(groupState->visible, groupState->selected, groupState->highlighted,
		groupState->lineStyle, groupState->lineWidth, groupState->transparency, groupState->colorOverride,
		groupState->color, groupState->materialColorValid, groupState->materialColor, groupState->materialRevision);
	    this->metadata->prepare();
	    if (this->metadata->changed()) this->hierarchy.effects().include(target, PeerEffects::Revision::Frame);
	}
	this->hierarchy.prepare();
	this->valid = true;
    }
    BObolSceneController &scene;
    std::string key, path;
    SbModernUtils::SoNodeRef sourceOwner;
    SoBRLDatabaseSource *source = nullptr;
    HierarchyPublication hierarchy;
    std::unique_ptr<SourceIndexPublication> indexes;
    std::unique_ptr<GroupPublication> metadata;
    bool valid = false, structural = false, published = false;
};

int
BObolSceneController::publishDatabaseSourceInstance(const BObolDatabaseSourcePublishState &state)
{
    return this->publishDatabaseSourceInstance(state, nullptr);
}

int
BObolSceneController::publishDatabaseSourceInstance(const BObolDatabaseSourcePublishState &state,
    const BObolSceneGroupPublishState *groupState)
{
    return this->publishDatabaseSourceInstanceImpl(nullptr, FALSE, state,
	groupState);
}

int
BObolSceneController::publishDatabaseSourceInstance(
    const char *currentSourceInstanceKey,
    const BObolDatabaseSourcePublishState &state,
    const BObolSceneGroupPublishState *groupState)
{
    if (!currentSourceInstanceKey || !currentSourceInstanceKey[0])
	return -1;
    return this->publishDatabaseSourceInstanceImpl(currentSourceInstanceKey,
	TRUE, state, groupState);
}

int
BObolSceneController::publishDatabaseSourceInstanceImpl(
    const char *currentSourceInstanceKey,
    SbBool requireExisting,
    const BObolDatabaseSourcePublishState &state,
    const BObolSceneGroupPublishState *groupState)
{
    BObolPerformanceTimer timer(BOBOL_PERF_SOURCE_REPLACE_US);
    if (timer.active()) bobol_performance_counter_add(BOBOL_PERF_SOURCE_REPLACE_CALLS, 1);
    if (!state.sourceInstanceKey || !state.sourceInstanceKey[0] || !state.sourcePath || !state.sourcePath[0])
	return -1;
    if (groupState && (!state.database || !state.targetGroupPath || !state.targetGroupPath[0])) return -1;
    if (!state.database)
	return this->removeDatabaseSourceInstance(requireExisting ?
	    currentSourceInstanceKey : state.sourceInstanceKey);
    if (!scene_root_group(this->d->root)) return -1;

    SourcePublication publication(*this, requireExisting ?
	currentSourceInstanceKey : state.sourceInstanceKey,
	requireExisting, state, groupState);
    if (!publication.isValid()) return -1;
    auto &source = publication.getSource();
    const uint32_t revision = state.sourceRevisionValid && state.sourceRevision ? state.sourceRevision :
	bobol_identity_successor_or_terminate(source.sourceRevision.getValue());
    int changed = 0;
    std::exception_ptr failure;
    try {
	changed = source.publishState(state, revision, SourcePublication::committed, &publication);
	if (!publication.didCommit() && publication.changesScene()) {
	    SourcePublication::committed(&publication);
	    changed = 1;
	}
    } catch (...) { failure = std::current_exception(); }
    publication.notify(failure);
    if (failure) std::rethrow_exception(failure);
    return changed;
}

int
BObolSceneController::subsumeDatabaseSourceInstances(
    const char *targetSourceInstanceKey,
    const char *const *sourceInstanceKeys,
    size_t sourceInstanceCount)
{
    if (!targetSourceInstanceKey || !targetSourceInstanceKey[0] ||
	(!sourceInstanceKeys && sourceInstanceCount) ||
	sourceInstanceCount > size_t(std::numeric_limits<int>::max()) ||
	!scene_root_group(this->d->root))
	return -1;

    SoBRLDatabaseSource *target =
	this->findIndexedDatabaseSourceInstance(targetSourceInstanceKey);
    if (!target)
	return 0;

    SceneChildOrders orders;
    std::vector<const SoBRLDatabaseSource *> sources;
    sources.reserve(sourceInstanceCount);
    std::unordered_set<SoBRLDatabaseSource *> selected;
    for (size_t i = 0; i < sourceInstanceCount; ++i) {
	const char *key = sourceInstanceKeys[i];
	if (!key || !key[0])
	    continue;
	SoBRLDatabaseSource *source =
	    this->findIndexedDatabaseSourceInstance(key);
	if (!source || source == target || !selected.insert(source).second)
	    continue;

	SoGroup *parent = this->findIndexedDatabaseSourceInstanceParent(key);
	int childIndex = parent ? scene_group_find_child_index(parent, source) : -1;
	if (!parent || childIndex < 0) {
	    parent = NULL;
	    childIndex = -1;
	    (void)find_database_source_instance_recursive(
		static_cast<SoGroup *>(this->d->root), key, &parent,
		&childIndex);
	}
	if (!parent || childIndex < 0)
	    continue;

	auto inserted = orders.emplace(parent, std::vector<SoNode *>());
	auto &children = inserted.first->second;
	if (inserted.second) {
	    children.reserve(size_t(parent->getNumChildren()));
	    for (int child = 0; child < parent->getNumChildren(); ++child)
		children.push_back(parent->getChild(child));
	}
	auto occurrence = std::find(children.begin(), children.end(), source);
	if (occurrence == children.end())
	    continue;
	children.erase(occurrence);
	sources.push_back(source);
    }
    if (sources.empty())
	return 0;

    class AdoptionEffects : public BObolSourceAdoptionEffects {
    public:
	AdoptionEffects(BObolSceneController &owner,
	    const SceneChildOrders &nextOrders, SoBRLDatabaseSource &destination) :
	    scene(owner), orders(nextOrders), target(destination)
	{}
	void prepare() override
	{
	    this->publication =
		std::make_unique<ChildPublication>(this->scene, this->orders, true);
	    this->publication->includeEffect(&this->target);
	}
	void commit() noexcept override { this->publication->commitPrepared(); }
	void notify() override
	{
	    std::exception_ptr failure;
	    this->publication->notify(failure);
	    if (failure) std::rethrow_exception(failure);
	}
    private:
	BObolSceneController &scene;
	const SceneChildOrders &orders;
	SoBRLDatabaseSource &target;
	std::unique_ptr<ChildPublication> publication;
    } effects(*this, orders, *target);

    (void)target->adoptCompactOccurrencesFrom(
	sources.data(), sources.size(), &effects);
    return static_cast<int>(sources.size());
}

namespace {

class SceneRealizationCommitWitness : public BObolSourceRealizationEffects {
public:
    explicit SceneRealizationCommitWitness(SbBool *result) : committed(result) {}
    void prepare(const SoBRLDatabaseSource &, bool, const SbString &,
	const std::vector<SoNode *> &, const std::vector<SoNode *> &) override {}
    void commit(bool changed) noexcept override
    {
	if (changed && this->committed)
	    *this->committed = TRUE;
    }
    void notify() override {}

private:
    SbBool *committed;
};

}

int
BObolSceneController::adoptDatabaseSourceInstanceRealization(
    const char *sourceInstanceKey,
    const BObolSourceRealizationStamp &stamp,
    SoBRLDatabaseSource *detached,
    SbBool authoritativeStreamDrained,
    const std::shared_ptr<BObolCompactOccurrenceStream> &stagedSourceStream,
    SbBool *publicationCommitted)
{
    if (publicationCommitted)
	*publicationCommitted = FALSE;
    if (!sourceInstanceKey || !sourceInstanceKey[0] || !detached ||
	!scene_root_group(this->d->root))
	return -1;

    SoBRLDatabaseSource *source =
	this->findIndexedDatabaseSourceInstance(sourceInstanceKey);
    if (!source || !source->matchesRealizationStamp(stamp))
	return 0;

    SbModernUtils::SoNodeRef sourceOwner(source);
    SbModernUtils::SoNodeRef detachedOwner(detached);
    SceneRealizationCommitWitness witness(publicationCommitted);
    auto effects = BObolSceneController::prepareRealizationEffects(*this,
	&witness);
    return source->adoptDetachedCompactRealization(stamp, detached,
	authoritativeStreamDrained, stagedSourceStream, effects.get());
}

int
BObolSceneController::mergeDatabaseSourceInstanceCompactOccurrences(
    const char *sourceInstanceKey,
    const BObolSourceRealizationStamp &stamp,
    const std::vector<BObolCompactOccurrence> &occurrences,
    SbBool authoritativeGeometry,
    SbBool *batchCompleted,
    size_t reserveCapacity)
{
    if (batchCompleted)
	*batchCompleted = FALSE;
    if (!sourceInstanceKey || !sourceInstanceKey[0] ||
	!scene_root_group(this->d->root))
	return -1;
    if (occurrences.empty())
	return 0;

    SoBRLDatabaseSource *source =
	this->findIndexedDatabaseSourceInstance(sourceInstanceKey);
    if (!source || !source->matchesRealizationStamp(stamp))
	return 0;

    SbModernUtils::SoNodeRef sourceOwner(source);
    PeerEffects effects(*this, source);
    return source->mergeCompactOccurrences(occurrences,
	authoritativeGeometry, PeerEffects::frameCommitted, &effects,
	batchCompleted, reserveCapacity);
}

int
BObolSceneController::certifyDatabaseSourceInstanceCompactStream(
    const char *sourceInstanceKey,
    const BObolSourceRealizationStamp &stamp,
    size_t expectedCount,
    const BObolCompactSourceProfile *profile)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0] || !expectedCount ||
	(profile && !profile->isValid(expectedCount)) ||
	!scene_root_group(this->d->root))
	return -1;

    SoBRLDatabaseSource *source =
	this->findIndexedDatabaseSourceInstance(sourceInstanceKey);
    if (!source || !source->matchesRealizationStamp(stamp))
	return 0;

    SbModernUtils::SoNodeRef sourceOwner(source);
    return source->certifyCompactStream(expectedCount, profile);
}

int
BObolSceneController::publishDatabaseSourceInstanceCompactSnapshot(
    const char *sourceInstanceKey,
    const BObolSourceRealizationStamp &stamp,
    const std::vector<BObolCompactOccurrence> &occurrences,
    const SbBox3f *certifiedBounds,
    const BObolCompactSourceProfile *profile,
    SbBool *publicationCommitted)
{
    if (publicationCommitted)
	*publicationCommitted = FALSE;
    if (!sourceInstanceKey || !sourceInstanceKey[0] || occurrences.empty() ||
	(certifiedBounds && certifiedBounds->isEmpty()) ||
	(profile && (!profile->isValid() ||
	 profile->occurrenceCount > std::numeric_limits<size_t>::max())) ||
	!scene_root_group(this->d->root))
	return -1;

    SoBRLDatabaseSource *source =
	this->findIndexedDatabaseSourceInstance(sourceInstanceKey);
    if (!source || !source->matchesRealizationStamp(stamp))
	return 0;

    SbModernUtils::SoNodeRef sourceOwner(source);
    SceneRealizationCommitWitness witness(publicationCommitted);
    auto effects = BObolSceneController::prepareRealizationEffects(*this,
	&witness);
    return source->publishCompactSnapshot(stamp, occurrences,
	certifiedBounds, profile, effects.get());
}

int
BObolSceneController::renameDatabaseSourceInstance(
    const char *sourceInstanceKey,
    const char *newSourceInstanceKey,
    const char *newSourcePath,
    uint32_t sourceRevision)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0] ||
	!newSourceInstanceKey || !newSourceInstanceKey[0] ||
	!newSourcePath || !newSourcePath[0])
	return -1;
    SoBRLDatabaseSource *existing = this->findDatabaseSourceInstance(sourceInstanceKey);
    if ((bu_strcmp(sourceInstanceKey, newSourceInstanceKey) == 0 ||
	 bu_strcmp(skip_leading_slash(sourceInstanceKey),
		skip_leading_slash(newSourceInstanceKey)) == 0) &&
	database_source_path_equal(existing, newSourcePath) &&
	(!sourceRevision || sourceRevision == existing->sourceRevision.getValue()))
	return 0;
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *sourceParent = NULL;
    int sourceIndex = -1;
    SoBRLDatabaseSource *source = existing;
    if (source) {
	sourceParent = this->findIndexedDatabaseSourceInstanceParent(
			   sourceInstanceKey);
	if (sourceParent)
	    sourceIndex = scene_group_find_child_index(sourceParent, source);
    }
    if (!source || !sourceParent || sourceIndex < 0) {
	SoGroup *rootGroup = static_cast<SoGroup *>(this->d->root);
	source = find_database_source_instance_recursive(rootGroup,
		 sourceInstanceKey, &sourceParent, &sourceIndex);
    }
    if (!source || !sourceParent || sourceIndex < 0)
	return 0;

    SoGroup *conflictParent = NULL;
    int conflictIndex = -1;
    SoBRLDatabaseSource *conflict =
	this->findIndexedDatabaseSourceInstance(newSourceInstanceKey);
    if (conflict) {
	conflictParent = this->findIndexedDatabaseSourceInstanceParent(
			     newSourceInstanceKey);
	if (conflictParent)
	    conflictIndex = scene_group_find_child_index(conflictParent,
			    conflict);
    }
    if (conflict == source) conflict = nullptr;
    if (conflict && (!conflictParent || conflictIndex < 0)) return -1;
    HierarchyPublication hierarchy(*this);
    SourceIndexPublication publication(*this, hierarchy, *source, *sourceParent, newSourceInstanceKey,
	newSourcePath, conflict, conflictParent, conflictIndex);
    if (!publication.isValid()) return -1;
    hierarchy.prepare();
    int changed = 0;
    std::exception_ptr failure;
    try {
	changed = source->retargetDatabaseSourceInstance(newSourceInstanceKey, newSourcePath,
	    sourceRevision, SourceIndexPublication::committed, &publication);
    } catch (...) { failure = std::current_exception(); }
    publication.notify(failure);
    if (failure) std::rethrow_exception(failure);
    return changed;
}

static bool
scene_retarget_source_instance_key(std::string &key,
    const char *currentSourcePath, const char *newSourcePath)
{
    const std::string current = normalized_index_key(currentSourcePath);
    const std::string next = normalized_index_key(newSourcePath);
    if (key.empty() || current.empty() || next.empty())
	return false;

    size_t position = key.rfind(current);
    while (position != std::string::npos) {
	const size_t end = position + current.size();
	const bool leftBoundary = position == 0 || key[position - 1] == '/' ||
	    key[position - 1] == ':';
	const bool rightBoundary = end == key.size() || key[end] == ':';
	if (leftBoundary && rightBoundary) {
	    key.replace(position, current.size(), next);
	    return true;
	}
	if (!position)
	    break;
	position = key.rfind(current, position - 1);
    }
    return false;
}

static SoGroup *
scene_source_parent_for_node(SoGroup *root, SoBRLDatabaseSource *target)
{
    if (!root || !target)
	return nullptr;
    std::vector<SoGroup *> pending{root};
    std::unordered_set<SoGroup *> visited;
    while (!pending.empty()) {
	SoGroup *parent = pending.back();
	pending.pop_back();
	if (!visited.insert(parent).second)
	    continue;
	for (int i = 0; i < parent->getNumChildren(); ++i) {
	    SoNode *child = parent->getChild(i);
	    if (child == target)
		return parent;
	    if (child && child->isOfType(SoGroup::getClassTypeId()))
		pending.push_back(static_cast<SoGroup *>(child));
	}
    }
    return nullptr;
}

static SbString
scene_retarget_group_intent_path(const SoBRLSceneGroup &group,
    const std::string &newGroupPath)
{
    const char *intentValue = group.drawIntentPath.getValue().getString();
    const std::string intent = intentValue ? intentValue : "";
    const std::string current = normalized_index_key(
	group.groupPath.getValue().getString());
    if (intent.empty() || current.empty() || intent.size() < current.size())
	return intent.c_str();
    const size_t position = intent.size() - current.size();
    if (intent.compare(position, current.size(), current) != 0 ||
	(position && intent[position - 1] != ':' &&
	 intent[position - 1] != '/'))
	return intent.c_str();
    std::string nextIntent = intent.substr(0, position);
    nextIntent += newGroupPath;
    return nextIntent.c_str();
}

int
BObolSceneController::renameDatabaseObject(
    const char *oldObjectPath, const char *newObjectPath,
    const std::vector<SbString> &pathDerivedSourceInstanceKeys,
    uint32_t sourceRevision)
{
    if (!oldObjectPath || !oldObjectPath[0] ||
	!newObjectPath || !newObjectPath[0] ||
	!scene_root_group(this->d->root))
	return -1;

    std::string identityProbe = normalized_index_key(oldObjectPath);
    if (!bobol_database_retarget_path_components(identityProbe,
	oldObjectPath, newObjectPath))
	return 0;

    if (this->d->activeRealizationAction)
	this->d->activeRealizationAction->stopSceneTraversal();
    if (!this->d->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();

    std::unordered_set<std::string> pathDerivedKeys;
    pathDerivedKeys.reserve(pathDerivedSourceInstanceKeys.size());
    for (const SbString &key : pathDerivedSourceInstanceKeys) {
	const std::string normalized = normalized_index_key(key.getString());
	if (normalized.empty() || !pathDerivedKeys.insert(normalized).second ||
	    !this->findIndexedDatabaseSourceInstance(normalized.c_str()))
	    return -1;
    }

    struct SourceChange {
	SoBRLDatabaseSource *source;
	SoGroup *parent;
	std::string currentKey;
	std::string nextKey;
	std::string nextPath;
	std::string nextParentKey;
	bool identityChanged;
	std::unique_ptr<BObolSourceRenamePublication> publication;
    };
    std::vector<SourceChange> plannedSourceChanges;
    plannedSourceChanges.reserve(
	this->d->indexes.databaseSourceOrder.size());
    std::unordered_set<std::string> consumedDerivedKeys;

    SoGroup *root = static_cast<SoGroup *>(this->d->root);
    for (SoBRLDatabaseSource *source :
	 this->d->indexes.databaseSourceOrder) {
	if (!source)
	    return -1;
	const SbString effectiveKey =
	    database_source_effective_instance_key(source);
	const std::string currentKey = effectiveKey.getString();
	std::string nextPath = source->path.getValue().getString();
	const bool pathChanged = bobol_database_retarget_path_components(
	    nextPath, oldObjectPath, newObjectPath);
	std::string nextKey = currentKey;
	const std::string normalizedKey = normalized_index_key(
	    currentKey.c_str());
	if (pathChanged && pathDerivedKeys.count(normalizedKey)) {
	    if (!scene_retarget_source_instance_key(nextKey,
		    source->path.getValue().getString(), nextPath.c_str()))
		return -1;
	    consumedDerivedKeys.insert(normalizedKey);
	}
	SoGroup *parent = this->findIndexedDatabaseSourceInstanceParent(
	    currentKey.c_str());
	if (!parent || scene_group_find_child_index(parent, source) < 0)
	    parent = scene_source_parent_for_node(root, source);
	if (!parent)
	    return -1;

	plannedSourceChanges.push_back({source, parent, currentKey, nextKey,
	    nextPath, source->parentInstanceKey.getValue().getString(),
	    pathChanged, nullptr});
    }
    if (consumedDerivedKeys.size() != pathDerivedKeys.size())
	return -1;

    std::unordered_map<std::string, std::string> renamedSourceKeys;
    renamedSourceKeys.reserve(plannedSourceChanges.size());
    for (const SourceChange &change : plannedSourceChanges) {
	const std::string current = normalized_index_key(
	    change.currentKey.c_str());
	const std::string next = normalized_index_key(change.nextKey.c_str());
	if (current != next)
	    renamedSourceKeys.emplace(current, change.nextKey);
    }

    std::vector<SourceChange> sourceChanges;
    sourceChanges.reserve(plannedSourceChanges.size());
    std::unordered_map<SoBRLDatabaseSource *, size_t> sourceChangeIndex;
    sourceChangeIndex.reserve(plannedSourceChanges.size());
    for (SourceChange &change : plannedSourceChanges) {
	const std::string currentParent = normalized_index_key(
	    change.nextParentKey.c_str());
	const auto renamedParent = renamedSourceKeys.find(currentParent);
	const bool parentChanged = renamedParent != renamedSourceKeys.end();
	if (parentChanged)
	    change.nextParentKey = renamedParent->second;
	change.publication = std::make_unique<BObolSourceRenamePublication>(
	    *change.source, oldObjectPath, newObjectPath,
	    change.identityChanged ? change.nextKey.c_str() : nullptr,
	    change.identityChanged ? change.nextPath.c_str() : nullptr,
	    sourceRevision,
	    parentChanged ? change.nextParentKey.c_str() : nullptr);
	if (!change.publication->changed())
	    continue;
	sourceChangeIndex.emplace(change.source, sourceChanges.size());
	sourceChanges.push_back(std::move(change));
    }

    struct GroupChange {
	SoBRLSceneGroup *group;
	std::string nextPath;
	SbString nextIntentPath;
	SbName oldName;
	SbName nextName;
	bool nameChanged;
	std::unique_ptr<GroupPublication> publication;
    };
    std::vector<GroupChange> groupChanges;
    std::unordered_map<SoBRLSceneGroup *, std::string> nextGroupPaths;
    std::unordered_set<SoBRLSceneGroup *> visitedGroups;
    groupChanges.reserve(this->d->indexes.groupPathIndex.size());
    nextGroupPaths.reserve(this->d->indexes.groupPathIndex.size());
    visitedGroups.reserve(this->d->indexes.groupPathIndex.size());
    for (const auto &entry : this->d->indexes.groupPathIndex) {
	SoGroup *indexed = entry.second;
	if (!indexed ||
	    !indexed->isOfType(SoBRLSceneGroup::getClassTypeId()))
	    continue;
	auto *group = static_cast<SoBRLSceneGroup *>(indexed);
	if (!visitedGroups.insert(group).second)
	    continue;
	std::string nextPath = group->groupPath.getValue().getString();
	if (!bobol_database_retarget_path_components(nextPath,
		oldObjectPath, newObjectPath))
	    continue;
	const std::string currentPath = normalized_index_key(
	    group->groupPath.getValue().getString());
	const size_t currentSlash = currentPath.find_last_of('/');
	const std::string currentLeaf = currentPath.substr(
	    currentSlash == std::string::npos ? 0 : currentSlash + 1);
	const size_t nextSlash = nextPath.find_last_of('/');
	const std::string nextLeaf = nextPath.substr(
	    nextSlash == std::string::npos ? 0 : nextSlash + 1);
	const bool nameChanged = currentLeaf != nextLeaf;
	GroupChange change{group, nextPath,
	    scene_retarget_group_intent_path(*group, nextPath),
	    group->getName(),
	    SbName(nameChanged ? nextLeaf.c_str() :
		group->getName().getString()), nameChanged, nullptr};
	change.publication = std::make_unique<GroupPublication>(*this,
	    *group);
	change.publication->setPath(change.nextPath.c_str());
	change.publication->setIntentPath(change.nextIntentPath);
	change.publication->prepare();
	nextGroupPaths.emplace(group, change.nextPath);
	groupChanges.push_back(std::move(change));
    }

    if (sourceChanges.empty() && groupChanges.empty())
	return 0;
    const size_t maximumChanges = static_cast<size_t>(
	std::numeric_limits<int>::max());
    if (sourceChanges.size() > maximumChanges ||
	groupChanges.size() > maximumChanges - sourceChanges.size())
	return -1;

    struct NameReservation {
	std::string name;
	SbModernUtils::SoNodeRef node;
    };
    std::vector<NameReservation> nameReservations;
    nameReservations.reserve(groupChanges.size() * 2);
    const auto reserveName = [&nameReservations](const SbName &value) {
	const std::string name = value.getString();
	if (std::any_of(nameReservations.begin(), nameReservations.end(),
		[&name](const NameReservation &reservation) {
		    return reservation.name == name;
		}))
	    return;
	auto *sentinel = new SoSeparator;
	SbModernUtils::SoNodeRef retainedSentinel(sentinel);
	sentinel->setName(value);
	nameReservations.push_back({name, std::move(retainedSentinel)});
    };
    for (const GroupChange &change : groupChanges) {
	if (!change.nameChanged)
	    continue;
	reserveName(change.oldName);
	reserveName(change.nextName);
    }
    std::vector<GroupChange *> preparedNames;
    preparedNames.reserve(groupChanges.size());
    try {
	for (GroupChange &change : groupChanges) {
	    if (!change.nameChanged)
		continue;
	    change.group->setName(change.nextName);
	    preparedNames.push_back(&change);
	}
    } catch (...) {
	for (auto current = preparedNames.rbegin();
	     current != preparedNames.rend(); ++current)
	    (*current)->group->setName((*current)->oldName);
	throw;
    }
    for (auto current = preparedNames.rbegin();
	 current != preparedNames.rend(); ++current)
	(*current)->group->setName((*current)->oldName);

    SceneIndexes nextIndexes;
    nextIndexes.groupPathIndex.reserve(
	this->d->indexes.groupPathIndex.size());
    visitedGroups.clear();
    for (const auto &entry : this->d->indexes.groupPathIndex) {
	SoGroup *group = entry.second;
	if (!group ||
	    !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	    continue;
	auto *sceneGroup = static_cast<SoBRLSceneGroup *>(group);
	if (!visitedGroups.insert(sceneGroup).second)
	    continue;
	const auto changed = nextGroupPaths.find(sceneGroup);
	const std::string path = changed == nextGroupPaths.end() ?
	    normalized_index_key(sceneGroup->groupPath.getValue().getString()) :
	    normalized_index_key(changed->second.c_str());
	if (path.empty() || !nextIndexes.groupPathIndex.emplace(path,
		group).second)
	    return -1;
    }

    const auto &order = this->d->indexes.databaseSourceOrder;
    nextIndexes.databaseSourceOrder = order;
    nextIndexes.databaseSourceOrderIndex.reserve(order.size());
    nextIndexes.databaseSourceRoutingIndex.reserve(order.size());
    nextIndexes.databaseSourceInstanceIndex.reserve(order.size());
    nextIndexes.databaseSourceInstanceParentIndex.reserve(order.size());
    nextIndexes.databaseSourcePathIndex.reserve(order.size());
    nextIndexes.databaseSourcePathInstancesIndex.reserve(order.size());
    for (size_t i = 0; i < order.size(); ++i) {
	SoBRLDatabaseSource *source = order[i];
	const auto changed = sourceChangeIndex.find(source);
	const SourceChange *change = changed == sourceChangeIndex.end() ?
	    nullptr : &sourceChanges[changed->second];
	const std::string key = normalized_index_key(change ?
	    change->nextKey.c_str() :
	    database_source_effective_instance_key(source).getString());
	const std::string path = normalized_index_key(change ?
	    change->nextPath.c_str() : source->path.getValue().getString());
	SoGroup *parent = change ? change->parent :
	    this->findIndexedDatabaseSourceInstanceParent(key.c_str());
	if (!parent)
	    parent = scene_source_parent_for_node(root, source);
	if (key.empty() || path.empty() || !parent ||
	    !nextIndexes.databaseSourceInstanceIndex.emplace(key,
		source).second ||
	    !nextIndexes.databaseSourceInstanceParentIndex.emplace(key,
		parent).second)
	    return -1;
	nextIndexes.databaseSourcePathIndex[path] = source;
	nextIndexes.databaseSourcePathInstancesIndex[path].push_back(source);
	if (!nextIndexes.databaseSourceRoutingIndex.emplace(
		source->getCompactSourceRoutingId(), source).second)
	    return -1;
	nextIndexes.databaseSourceOrderIndex.emplace(source, i);
    }

    struct RepositoryRename {
	std::shared_ptr<BObolRealizationRepository> repository;
	std::unique_ptr<BObolRealizationRepository::ObjectRename> publication;
    };
    std::vector<RepositoryRename> repositoryRenames;
    std::unordered_set<BObolRealizationRepository *> preparedRepositories;
    std::set<BObolSceneController *> owners;
    const auto includeOwners = [&owners](SoNode *node) {
	RootOwnership::visitControllers(node, nullptr,
	    [&owners](BObolSceneController *owner) {
		owners.insert(owner);
	    });
    };
    for (const SourceChange &change : sourceChanges)
	includeOwners(change.source);
    for (const GroupChange &change : groupChanges)
	includeOwners(change.group);
    owners.insert(this);
    repositoryRenames.reserve(owners.size());
    preparedRepositories.reserve(owners.size());
    for (BObolSceneController *owner : owners) {
	if (owner->d->activeRealizationAction)
	    owner->d->activeRealizationAction->stopSceneTraversal();
	auto repository = owner->d->realizationRepository;
	if (!repository ||
	    !preparedRepositories.insert(repository.get()).second)
	    continue;
	auto publication = repository->prepareObjectRename(oldObjectPath,
	    newObjectPath);
	if (publication)
	    repositoryRenames.push_back({std::move(repository),
		std::move(publication)});
    }
    PeerEffects effects(*this);
    bool structural = !groupChanges.empty();
    for (const SourceChange &change : sourceChanges) {
	structural = structural || change.identityChanged;
	effects.include(change.source, change.identityChanged ?
	    PeerEffects::Revision::Structural : PeerEffects::Revision::Frame);
    }
    for (const GroupChange &change : groupChanges)
	effects.include(change.group, PeerEffects::Revision::Structural);

    /* The retained sentinels keep the preflighted Coin name-list capacity
     * alive, so these non-notifying assignments cannot allocate. */
    for (GroupChange &change : groupChanges) {
	if (change.nameChanged)
	    change.group->setName(change.nextName);
    }
    for (SourceChange &change : sourceChanges)
	change.publication->commit();
    for (GroupChange &change : groupChanges)
	change.publication->commit();
    for (RepositoryRename &rename : repositoryRenames)
	rename.publication->commit();
    this->d->indexes.swap(nextIndexes);
    this->d->databaseSourceIndexValid = TRUE;
    effects.commit();
    if (structural)
	this->advanceStructuralRevision();
    else
	this->advanceFrameRevision();

    for (SourceChange &change : sourceChanges)
	change.publication->restore();
    for (GroupChange &change : groupChanges)
	change.publication->restore();
    std::exception_ptr failure;
    for (SourceChange &change : sourceChanges)
	change.publication->notify(failure);
    for (GroupChange &change : groupChanges)
	change.publication->notify(failure);
    if (failure)
	std::rethrow_exception(failure);
    return static_cast<int>(sourceChanges.size() + groupChanges.size());
}

int
BObolSceneController::setDatabaseSourceState(const char *sourcePath,
	SbBool sourceRevisionValid,
	uint32_t sourceRevision,
	uint32_t inputsRevision,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceState(sourceInstanceKey.getString(),
	    sourceRevisionValid, sourceRevision, inputsRevision, visible,
	    selected, highlighted, lineStyle, lineWidth, transparency, colorOverride,
	    color, materialColorValid, materialColor, materialRevision);
}

int
BObolSceneController::setDatabaseSourceInstanceState(
    const char *sourceInstanceKey,
    SbBool sourceRevisionValid,
    uint32_t sourceRevision,
    uint32_t inputsRevision,
    SbBool visible,
    SbBool selected,
    SbBool highlighted,
    int lineStyle,
    int lineWidth,
    float transparency,
    SbBool colorOverride,
    const SbColor &color,
    SbBool materialColorValid,
    const SbColor &materialColor,
    uint32_t materialRevision)
{
    BObolPerformanceTimer timer(BOBOL_PERF_SOURCE_STATE_US);
    if (timer.active())
	bobol_performance_counter_add(BOBOL_PERF_SOURCE_STATE_CALLS, 1);

    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    PeerEffects effects(*this, source);
    return source->setDisplayState(sourceRevisionValid,
			sourceRevision, inputsRevision, visible, selected, highlighted, lineStyle,
			lineWidth, transparency, colorOverride, color, materialColorValid,
			materialColor, materialRevision, PeerEffects::frameCommitted, &effects);
}

int
BObolSceneController::setDatabaseSourceDisplayPatch(const char *sourcePath,
	const BObolDatabaseSourceDisplayPatch &patch)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceDisplayPatch(
	       sourceInstanceKey.getString(), patch);
}

int
BObolSceneController::setDatabaseSourceInstanceDisplayPatch(
    const char *sourceInstanceKey,
    const BObolDatabaseSourceDisplayPatch &patch)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    PeerEffects effects(*this, source);
    return source->applyDisplayPatch(patch, PeerEffects::frameCommitted, &effects);
}

static bool
scene_group_display_patch_supported(
    const BObolDatabaseSourceDisplayPatch &patch)
{
    return !patch.selectedColorValid && !patch.highlightedColorValid &&
	!patch.ghostedColorValid;
}

static BObolDatabaseSourceDisplayPatch
scene_group_display_patch_snapshot(const SoBRLSceneGroup &group,
    const BObolDatabaseSourceDisplayPatch &requested)
{
    BObolDatabaseSourceDisplayPatch snapshot = requested;
    if (snapshot.visibleValid) snapshot.visible = group.visible.getValue();
    if (snapshot.selectedValid) snapshot.selected = group.selected.getValue();
    if (snapshot.highlightedValid)
	snapshot.highlighted = group.highlighted.getValue();
    if (snapshot.lineStyleValid) snapshot.lineStyle = group.lineStyle.getValue();
    if (snapshot.lineWidthValid) snapshot.lineWidth = group.lineWidth.getValue();
    if (snapshot.transparencyValid)
	snapshot.transparency = group.transparency.getValue();
    if (snapshot.colorOverrideValid)
	snapshot.colorOverride = group.colorOverride.getValue();
    if (snapshot.colorValid) snapshot.color = group.color.getValue();
    return snapshot;
}

static bool
scene_group_display_patch_matches(const SoBRLSceneGroup &group,
    const BObolDatabaseSourceDisplayPatch &snapshot)
{
    return (!snapshot.visibleValid ||
	group.visible.getValue() == snapshot.visible) &&
	(!snapshot.selectedValid ||
	 group.selected.getValue() == snapshot.selected) &&
	(!snapshot.highlightedValid ||
	 group.highlighted.getValue() == snapshot.highlighted) &&
	(!snapshot.lineStyleValid ||
	 group.lineStyle.getValue() == snapshot.lineStyle) &&
	(!snapshot.lineWidthValid ||
	 group.lineWidth.getValue() == snapshot.lineWidth) &&
	(!snapshot.transparencyValid ||
	 !scene_group_float_different(group.transparency.getValue(),
	     snapshot.transparency)) &&
	(!snapshot.colorOverrideValid ||
	 group.colorOverride.getValue() == snapshot.colorOverride) &&
	(!snapshot.colorValid ||
	 scene_group_color_equal(group.color.getValue(), snapshot.color));
}

int
BObolSceneController::applyPresentationTransaction(
    const BObolScenePresentationTransaction &transaction)
{
    struct SourceTarget {
	SbModernUtils::SoNodeRef owner;
	SbString key;
	BObolDatabaseSourcePresentationPatch patch;
	uint64_t cadRevision;
	SbUniqueId nodeId;
    };
    struct GroupTarget {
	SbModernUtils::SoNodeRef owner;
	SbString path;
	BObolDatabaseSourceDisplayPatch requested;
	BObolDatabaseSourceDisplayPatch snapshot;
    };

    std::vector<SourceTarget> sources;
    std::vector<GroupTarget> groups;
    std::set<std::string> sourceKeys;
    std::set<std::string> groupPaths;
    sources.reserve(transaction.sources.size());
    groups.reserve(transaction.groups.size());

    for (const auto &requested : transaction.sources) {
	const std::string key = normalized_index_key(
	    requested.sourceInstanceKey.getString());
	if (key.empty() || !sourceKeys.insert(key).second)
	    return -1;
	if (!SoBRLDatabaseSource::presentationPatchValid(
		requested.presentation))
	    return -1;
	auto *source = this->findDatabaseSourceInstance(key.c_str());
	if (!source) {
	    if (requested.expectedStamp.isValid())
		continue;
	    return -1;
	}
	if (requested.expectedStamp.isValid() &&
		!source->matchesPresentationStamp(requested.expectedStamp))
	    continue;
	sources.push_back({SbModernUtils::SoNodeRef(source), key.c_str(),
	    requested.presentation, source->cadBatchRevisionGet(),
	    source->getNodeId()});
    }

    for (const auto &requested : transaction.groups) {
	const std::string path = normalized_index_key(
	    requested.groupPath.getString());
	if (path.empty() || !groupPaths.insert(path).second ||
	    !scene_group_display_patch_supported(requested.display))
	    return -1;
	SoGroup *node = this->findGroup(path.c_str());
	if (!node || !node->isOfType(SoBRLSceneGroup::getClassTypeId()))
	    return -1;
	auto *group = static_cast<SoBRLSceneGroup *>(node);
	groups.push_back({SbModernUtils::SoNodeRef(group), path.c_str(),
	    requested.display,
	    scene_group_display_patch_snapshot(*group, requested.display)});
    }

    int changed = 0;
    for (auto &target : sources) {
	auto *source = static_cast<SoBRLDatabaseSource *>(target.owner.get());
	if (this->findDatabaseSourceInstance(target.key.getString()) != source ||
	    source->cadBatchRevisionGet() != target.cadRevision ||
	    source->getNodeId() != target.nodeId)
	    continue;
	PeerEffects effects(*this, source);
	const int result = source->applyPresentationPatch(target.patch,
	    PeerEffects::frameCommitted, &effects);
	if (result < 0)
	    throw std::logic_error(
		"validated source presentation patch was rejected");
	if (result > 0)
	    ++changed;
    }

    for (auto &target : groups) {
	auto *group = static_cast<SoBRLSceneGroup *>(target.owner.get());
	if (this->findGroup(target.path.getString()) != group ||
	    !scene_group_display_patch_matches(*group, target.snapshot))
	    continue;
	GroupPublication publication(*this, *group);
	publication.setDisplayPatch(target.requested);
	if (publication.publish() > 0)
	    ++changed;
    }
    return changed;
}

int
BObolSceneController::setDatabaseSourceDisplayName(const char *sourcePath,
	const char *displayName)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceDisplayName(
	       sourceInstanceKey.getString(),
	       displayName);
}

int
BObolSceneController::setDatabaseSourceInstanceDisplayName(
    const char *sourceInstanceKey,
    const char *displayName)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    PeerEffects effects(*this, source);
    return source->setDisplayNameState(displayName, PeerEffects::frameCommitted, &effects);
}

int
BObolSceneController::setDatabaseSourceDrawMode(const char *sourcePath,
	int drawMode)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceDrawMode(
	       sourceInstanceKey.getString(), drawMode);
}

int
BObolSceneController::setDatabaseSourcesEvaluatedRegionForPath(SoBRLDatabaseSource *const *sources,
    size_t count, const char *path, SbBool evaluated)
{
    if (!count)
	return 0;
    if (!sources)
	return -1;
    const SbString query(path ? path : "");
    struct Target {
	SbModernUtils::SoNodeRef owner;
	SbString key;
	uint64_t revision;
	SbUniqueId nodeId;
    };
    std::vector<Target> targets;
    targets.reserve(count);
    for (size_t i = 0; i < count; ++i) {
	auto *source = sources[i];
	const SbString key = database_source_effective_instance_key(source);
	if (!source || this->findDatabaseSourceInstance(key.getString()) != source)
	    return -1;
	targets.push_back({SbModernUtils::SoNodeRef(source), key,
	    source->cadBatchRevisionGet(), source->getNodeId()});
    }
    int changed = 0;
    for (auto &target : targets) {
	auto *source = static_cast<SoBRLDatabaseSource *>(target.owner.get());
	if (this->findDatabaseSourceInstance(target.key.getString()) != source ||
	    source->cadBatchRevisionGet() != target.revision || source->getNodeId() != target.nodeId)
	    continue;
	PeerEffects effects(*this, source);
	if (source->setEvaluatedRegionForPath(query.getString(), evaluated, PeerEffects::frameCommitted, &effects) > 0)
	    ++changed;
    }
    return changed;
}

int
BObolSceneController::setDatabaseSourceInstanceDrawMode(const char *sourceInstanceKey, int drawMode)
{
    SoBRLDatabaseSource *source = this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source) return -1;

    const int mode = drawMode == SoBRLDatabaseSource::SHADED ? SoBRLDatabaseSource::SHADED : SoBRLDatabaseSource::WIREFRAME;
    int representation = source->representationMode.getValue();
    const char *key = source->representationKey.getValue().getString();
    if (representation == SoBRLDatabaseSource::REPRESENTATION_DEFAULT ||
	representation == SoBRLDatabaseSource::REPRESENTATION_WIRE ||
	representation == SoBRLDatabaseSource::REPRESENTATION_SHADED) {
	representation = mode == SoBRLDatabaseSource::SHADED ?
	    SoBRLDatabaseSource::REPRESENTATION_SHADED : SoBRLDatabaseSource::REPRESENTATION_WIRE;
	if (!key[0]) key = source->instanceKey.getValue().getLength() ?
	    source->instanceKey.getValue().getString() : source->path.getValue().getString();
    }
    PeerEffects effects(*this, source);
    return source->publishDrawRepresentationState(mode, key, representation,
	SoBRLDatabaseSource::ConfigurationIntent::Draw, PeerEffects::frameCommitted, &effects);
}

int
BObolSceneController::setDatabaseSourceInstanceRepresentation(const char *sourceInstanceKey,
    const char *sourceRepresentationKey, int sourceRepresentationMode)
{
    SoBRLDatabaseSource *source = this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source) return -1;
    PeerEffects effects(*this, source);
    return source->publishDrawRepresentationState(source->drawMode.getValue(), sourceRepresentationKey,
	sourceRepresentationMode, SoBRLDatabaseSource::ConfigurationIntent::Representation, PeerEffects::frameCommitted, &effects);
}

int
BObolSceneController::setDatabaseSourceMaterialPolicy(const char *sourcePath,
	int materialPolicy)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceMaterialPolicy(
	       sourceInstanceKey.getString(),
	       materialPolicy);
}

int
BObolSceneController::setDatabaseSourceInstanceMaterialPolicy(
    const char *sourceInstanceKey,
    int materialPolicy)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    PeerEffects effects(*this, source);
    return source->setMaterialPolicyState(materialPolicy, PeerEffects::frameCommitted, &effects);
}

int
BObolSceneController::applyDatabaseSourceInstanceDrawMetadata(const char *sourceInstanceKey,
    const BObolDrawMetadataRecord &metadata, const char *path)
{
    if (!metadata.directoryFound)
	return -1;
    auto *source = this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;
    constexpr float colorMaximum = 255.0f;
    const SbColor color = metadata.hasColor ? SbColor(metadata.color[0] / colorMaximum,
	metadata.color[1] / colorMaximum, metadata.color[2] / colorMaximum) : SbColor(1, 1, 1);
    const SbString shader(metadata.hasShader ? metadata.shader : "");
    const int region = metadata.hasRegionId ? metadata.regionId : 0;
    const int air = metadata.hasAircode ? metadata.aircode : 0;
    const int material = metadata.hasMaterialId ? metadata.materialId : 0;
    const int los = metadata.hasLos ? metadata.los : 0;
    PeerEffects effects(*this, source);
    if (path && source->hasCompactInstanceIndex())
	return source->setCompactInstanceMetadataForPath(path, TRUE, region, air, material, los,
	    metadata.hasColor ? TRUE : FALSE, color, shader, PeerEffects::frameCommitted, &effects);
    return source->setDatabaseMetadataState(TRUE, region, air, material, los,
	metadata.hasColor ? TRUE : FALSE, color, shader, PeerEffects::frameCommitted, &effects);
}

bool
BObolSceneController::databaseSourcePublicationAccepted(SoBRLDatabaseSource *source, void *context)
{
    auto &scene = *static_cast<BObolSceneController *>(context);
    return source && scene.findDatabaseSourceInstance(database_source_effective_instance_key(source).getString()) == source;
}

int
BObolSceneController::refreshDatabaseSourceInstanceMaterialColorFromDatabase(
    const char *sourceInstanceKey,
    uint32_t materialRevision,
    struct db_i *overrideDbip)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    PeerEffects effects(*this, source);
    return source->refreshMaterialColorFromDatabase(materialRevision, overrideDbip,
	PeerEffects::frameCommitted, &effects);
}

int
BObolSceneController::refreshDatabaseSourceMaterialColorsFromDatabase(
    uint32_t materialRevision,
    struct db_i *overrideDbip)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    const int sourceCount = this->getDatabaseSourceCount();
    std::vector<SoBRLDatabaseSource *> sources;
    sources.reserve(static_cast<size_t>(sourceCount));
    for (int i = 0; i < sourceCount; i++) {
	SoBRLDatabaseSource *source = this->getDatabaseSource(i);
	if (source)
	    sources.push_back(source);
    }
    auto *database = overrideDbip;
    if (!database && !sources.empty())
	database = sources.front()->getDatabase();
    if (!database && !sources.empty())
	return -1;
    return this->refreshDatabaseSourcesMaterialColorsFromDatabase(sources.data(), sources.size(),
	materialRevision, database);
}

int
BObolSceneController::refreshDatabaseSourcesMaterialColorsFromDatabase(SoBRLDatabaseSource *const *sources,
    size_t count, uint32_t materialRevision, struct db_i *overrideDbip)
{
    if (!count)
	return 0;
    if (!sources)
	return -1;
    for (size_t i = 0; i < count; ++i)
	if (!databaseSourcePublicationAccepted(sources[i], this))
	    return -1;
    struct TargetEffects {
	BObolSceneController &scene;
	std::unique_ptr<PeerEffects> current;
    } effects{*this, nullptr};
    const auto accepted = [](SoBRLDatabaseSource *source, void *context) {
	auto &state = *static_cast<TargetEffects *>(context);
	if (!databaseSourcePublicationAccepted(source, &state.scene)) return false;
	// The sweep accepts each target after the preceding target's callbacks.
	// Prepare only its current owners, preserving completed-prefix effects.
	state.current = std::make_unique<PeerEffects>(state.scene, source);
	return true;
    };
    const auto committed = [](void *context) noexcept {
	PeerEffects::frameCommitted(static_cast<TargetEffects *>(context)->current.get());
    };
    return SoBRLDatabaseSource::refreshMaterialColors(sources, count, materialRevision, overrideDbip,
	committed, &effects, accepted);
}

int
BObolSceneController::setDatabaseSourcePlacementState(const char *sourcePath,
	SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstancePlacementState(
	       sourceInstanceKey.getString(),
	       drawMatrixValid, drawMatrix, drawCenterValid, drawCenter,
	       drawSizeValid, drawSize);
}

int
BObolSceneController::setDatabaseSourceInstancePlacementState(
    const char *sourceInstanceKey,
    SbBool drawMatrixValid,
    const SbMatrix &drawMatrix,
    SbBool drawCenterValid,
    const SbVec3f &drawCenter,
    SbBool drawSizeValid,
    float drawSize)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    PeerEffects effects(*this, source);
    return source->setPlacementState(drawMatrixValid,
	drawMatrix, drawCenterValid, drawCenter, drawSizeValid, drawSize,
	PeerEffects::frameCommitted, &effects);
}

int
BObolSceneController::setDatabaseSourceInstanceHierarchyState(
    const char *sourceInstanceKey,
    const char *parentInstanceKey,
    uint32_t occurrenceIndex,
    int booleanOperation)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    PeerEffects effects(*this, source);
    return source->setHierarchyState(parentInstanceKey, occurrenceIndex, booleanOperation,
	PeerEffects::frameCommitted, &effects);
}

int
BObolSceneController::setDatabaseSourceBoundsState(const char *sourcePath,
	SbBool boundsValid,
	const SbVec3f &boundsMin,
	const SbVec3f &boundsMax,
	SbBool boundsExact)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceBoundsState(
	       sourceInstanceKey.getString(), boundsValid,
	       boundsMin, boundsMax, boundsExact);
}

int
BObolSceneController::setDatabaseSourceInstanceBoundsState(
    const char *sourceInstanceKey,
    SbBool boundsValid,
    const SbVec3f &boundsMin,
    const SbVec3f &boundsMax,
    SbBool boundsExact,
    SbBool *publicationCommitted)
{
    if (publicationCommitted)
	*publicationCommitted = FALSE;
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    PeerEffects effects(*this, source);
    struct CommitWitness {
	PeerEffects &effects;
	SbBool *committed;
	static void publish(void *context) noexcept
	{
	    auto &self = *static_cast<CommitWitness *>(context);
	    PeerEffects::frameCommitted(&self.effects);
	    if (self.committed)
		*self.committed = TRUE;
	}
    } witness{effects, publicationCommitted};
    return source->setSourceBoundsState(boundsValid, boundsMin, boundsMax, boundsExact,
	CommitWitness::publish, &witness);
}

int
BObolSceneController::markDatabaseSourceStale(const char *sourcePath,
	uint32_t staleReason)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->markDatabaseSourceInstanceStale(sourceInstanceKey.getString(),
	    staleReason);
}

int
BObolSceneController::markDatabaseSourceInstanceStale(
    const char *sourceInstanceKey,
    uint32_t staleReason)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    if (!staleReason)
	staleReason = SoBRLDatabaseSource::STALE_SOURCE;

    if (this->d->realizationRepository) {
	if (staleReason & (SoBRLDatabaseSource::STALE_DATABASE |
		SoBRLDatabaseSource::STALE_INPUTS |
		SoBRLDatabaseSource::STALE_TESSELLATION)) {
	    this->d->realizationRepository->clear();
	} else if (staleReason & SoBRLDatabaseSource::STALE_VIEW) {
	    this->d->realizationRepository->invalidateViewVariants();
	} else if (staleReason & SoBRLDatabaseSource::STALE_SOURCE) {
	    const char *path = source->path.getValue().getString();
	    const char *name = path ? strrchr(path, '/') : NULL;
	    name = name && name[1] ? name + 1 : path;
	    struct directory *dp = source->getDatabase() && name && name[0] ?
		db_lookup(source->getDatabase(), name, LOOKUP_QUIET) : NULL;
	    if (dp && !(dp->d_flags & RT_DIR_COMB))
		this->d->realizationRepository->invalidateObject(name);
	}
    }

    uint32_t nextSourceRevision = source->sourceRevision.getValue();
    if (staleReason & (SoBRLDatabaseSource::STALE_SOURCE |
		       SoBRLDatabaseSource::STALE_DATABASE))
	bobol_identity_advance(nextSourceRevision);

    const uint32_t nextReason =
	source->staleReason.getValue() | staleReason;
    const int changed =
	source->sourceRevision.getValue() != nextSourceRevision ||
	!source->stale.getValue() ||
	source->staleReason.getValue() != nextReason ||
	source->realizationStatus.getValue() != SoBRLDatabaseSource::UNREALIZED ||
	source->realizationDiagnostic.getValue().getLength() > 0;
    if (!changed)
	return 0;

    source->sourceRevision = nextSourceRevision;
    source->markStale(staleReason);
    this->advanceFrameRevision();
    return 1;
}

int
BObolSceneController::refreshDatabaseSourceInstanceObject(
    const char *sourceInstanceKey, const char *objectPath,
    uint32_t sourceRevision)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	source = this->findDatabaseSource(sourceInstanceKey);
    if (!source)
	return -1;
    if (this->d->realizationRepository && objectPath && objectPath[0]) {
	const char *name = strrchr(objectPath, '/');
	name = name && name[1] ? name + 1 : objectPath;
	struct directory *dp = source->getDatabase() && name[0] ?
	    db_lookup(source->getDatabase(), name, LOOKUP_QUIET) : NULL;
	if (!dp || !(dp->d_flags & RT_DIR_COMB))
	    this->d->realizationRepository->invalidateObject(name);
    }
    const int changed = source->refreshCompactObjectGeometry(objectPath,
	sourceRevision);
    if (changed > 0) {
	if (this->d->realizationRepository)
	    this->d->realizationRepository->seedSource(source);
	this->advanceFrameRevision();
    }
    return changed;
}

int
BObolSceneController::setDatabaseSourceRealizationState(const char *sourcePath,
	int realizationStatus,
	uint32_t realizedSourceRevision,
	uint32_t realizedInputsRevision,
	uint32_t staleReason,
	const char *diagnostic)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceRealizationState(
	       sourceInstanceKey.getString(),
	       realizationStatus, realizedSourceRevision, realizedInputsRevision,
	       staleReason, diagnostic);
}

int
BObolSceneController::setDatabaseSourceInstanceRealizationState(
    const char *sourceInstanceKey,
    int realizationStatus,
    uint32_t realizedSourceRevision,
    uint32_t realizedInputsRevision,
    uint32_t staleReason,
    const char *diagnostic)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    PeerEffects effects(*this, source);
    return source->publishRealizationState(realizationStatus,
	realizedSourceRevision, realizedInputsRevision, staleReason, diagnostic,
	nullptr, PeerEffects::frameCommitted, &effects);
}

int
BObolSceneController::setDatabaseSourceInstanceRealizationState(
    const char *sourceInstanceKey,
    int realizationStatus,
    uint32_t realizedSourceRevision,
    uint32_t realizedInputsRevision,
    uint32_t staleReason,
    const char *diagnostic,
    int roleFlags)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    PeerEffects effects(*this, source);
    return source->publishRealizationStateWithRoles(realizationStatus,
	realizedSourceRevision, realizedInputsRevision, staleReason, diagnostic,
	roleFlags, nullptr, PeerEffects::frameCommitted, &effects);
}

int
BObolSceneController::setDatabaseSourceRealizationRoleFlags(
    const char *sourcePath,
    int roleFlags)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceRealizationRoleFlags(
	       sourceInstanceKey.getString(),
	       roleFlags);
}

int
BObolSceneController::setDatabaseSourceInstanceRealizationRoleFlags(
    const char *sourceInstanceKey,
    int roleFlags)
{
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    PeerEffects effects(*this, source);
    return source->setRealizationRoleFlags(roleFlags, PeerEffects::frameCommitted, &effects);
}

int
BObolSceneController::setDatabaseSourceRealizationViewPolicy(
    const char *sourcePath,
    SbBool viewDependent,
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
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SbString sourceInstanceKey =
	this->databaseSourceInstanceKeyForPath(sourcePath);
    if (sourceInstanceKey.getLength() == 0)
	return -1;

    return this->setDatabaseSourceInstanceRealizationViewPolicy(
	       sourceInstanceKey.getString(),
	       viewDependent, csgLodEnabled, meshLodEnabled,
	       viewScale, lodScale, viewWidth, viewHeight, botThreshold,
	       curveScale, pointScale);
}

int
BObolSceneController::setDatabaseSourceInstanceRealizationViewPolicy(
    const char *sourceInstanceKey,
    SbBool viewDependent,
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
    SoBRLDatabaseSource *source =
	this->findDatabaseSourceInstance(sourceInstanceKey);
    if (!source)
	return -1;

    PeerEffects effects(*this, source);
    return source->setRealizationViewPolicy(viewDependent, csgLodEnabled, meshLodEnabled,
	viewScale, lodScale, viewWidth, viewHeight, botThreshold, curveScale, pointScale,
	PeerEffects::frameCommitted, &effects);
}

int
BObolSceneController::moveDatabaseSourceToGroup(const char *sourcePath,
	const char *groupPath)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoBRLDatabaseSource *source = this->findIndexedDatabaseSource(sourcePath);
    SbString sourceInstanceKey = database_source_effective_instance_key(source);
    if (sourceInstanceKey.getLength() == 0)
	return 0;

    return this->moveDatabaseSourceInstanceToGroup(
	       sourceInstanceKey.getString(), groupPath);
}

int
BObolSceneController::moveDatabaseSourceInstanceToGroup(
    const char *sourceInstanceKey,
    const char *groupPath)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0] || !groupPath)
	return -1;
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *sourceParent = NULL;
    int sourceIndex = -1;
    SoBRLDatabaseSource *source =
	this->findIndexedDatabaseSourceInstance(sourceInstanceKey);
    if (source) {
	sourceParent = this->findIndexedDatabaseSourceInstanceParent(
			   sourceInstanceKey);
	if (sourceParent)
	    sourceIndex = scene_group_find_child_index(sourceParent, source);
    }
    if (!source || !sourceParent || sourceIndex < 0) {
	SoGroup *rootGroup = static_cast<SoGroup *>(this->d->root);
	source = find_database_source_instance_recursive(rootGroup,
		 sourceInstanceKey, &sourceParent, &sourceIndex);
    }
    if (!source || !sourceParent || sourceIndex < 0)
	return 0;

    SbString sourceParentPath = scene_group_index_path(sourceParent);
    if (sourceParentPath.getLength() == 0 &&
	sourceParent == scene_root_group(this->d->root))
	sourceParentPath = "/";
    if (scene_path_equal(sourceParentPath.getString(), groupPath))
	return 0;

    BObolPerformanceTimer timer(BOBOL_PERF_SOURCE_MOVE_US);
    if (timer.active())
	bobol_performance_counter_add(BOBOL_PERF_SOURCE_MOVE_CALLS, 1);

    SourcePublication publication(*this, *source, *sourceParent, groupPath);
    if (!publication.isValid()) return -1;
    if (!publication.changesScene()) return 0;
    SourcePublication::committed(&publication);
    std::exception_ptr failure;
    publication.notify(failure);
    if (failure) std::rethrow_exception(failure);
    return 1;
}

int
BObolSceneController::removeDatabaseSource(const char *sourcePath)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoBRLDatabaseSource *source = this->findIndexedDatabaseSource(sourcePath);
    return ChildPublication::removeSource(*this, source);
}

int
BObolSceneController::removeDatabaseSourceInstance(
    const char *sourceInstanceKey)
{
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return 0;
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *sourceParent = NULL;
    int childIndex = -1;
    SoBRLDatabaseSource *source =
	this->findIndexedDatabaseSourceInstance(sourceInstanceKey);
    if (source) {
	sourceParent = this->findIndexedDatabaseSourceInstanceParent(
			   sourceInstanceKey);
	if (sourceParent)
	    childIndex = scene_group_find_child_index(sourceParent, source);
    }
    if (!source || !sourceParent || childIndex < 0) {
	SoGroup *group = static_cast<SoGroup *>(this->d->root);
	source = find_database_source_instance_recursive(group, sourceInstanceKey,
		 &sourceParent, &childIndex);
    }
    if (!sourceParent || childIndex < 0)
	return 0;

    ChildPublication publication(*this, *sourceParent, *source, false);
    return publication.publish();
}

int
BObolSceneController::applyRemovalTransaction(
    const BObolSceneRemovalTransaction &transaction)
{
    SoGroup *root = scene_root_group(this->d->root);
    if (!root)
	return -1;

    struct GroupTarget {
	std::string path;
	SoGroup *group;
	SoGroup *parent;
    };
    std::vector<GroupTarget> groups;
    std::set<std::string> requestedGroupPaths;
    for (const SbString &requested : transaction.groupPaths) {
	const std::string path = normalized_index_key(requested.getString());
	if (path.empty())
	    return -1;
	if (!requestedGroupPaths.insert(path).second)
	    continue;

	SoGroup *group = this->findIndexedGroup(path.c_str());
	if (!group || group == root)
	    continue;
	std::vector<std::string> components;
	scene_group_path_components(path.c_str(), components);
	SoGroup *parent = root;
	for (size_t i = 0; parent && i + 1 < components.size(); ++i)
	    parent = scene_group_find_child(parent, components[i].c_str());
	if (!parent || scene_group_find_child_index(parent, group) < 0)
	    continue;
	groups.push_back({path, group, parent});
    }

    std::sort(groups.begin(), groups.end(), [](const GroupTarget &left,
	const GroupTarget &right) {
	if (left.path.size() != right.path.size())
	    return left.path.size() < right.path.size();
	return left.path < right.path;
    });
    std::vector<GroupTarget> roots;
    roots.reserve(groups.size());
    for (const GroupTarget &candidate : groups) {
	const bool covered = std::any_of(roots.begin(), roots.end(),
	    [&candidate](const GroupTarget &selected) {
		return scene_path_contains(selected.path, candidate.path);
	    });
	if (!covered)
	    roots.push_back(candidate);
    }

    struct SourceTarget {
	SoBRLDatabaseSource *source;
	SoGroup *parent;
    };
    std::vector<SourceTarget> sources;
    std::set<std::string> requestedSourceKeys;
    std::unordered_set<SoBRLDatabaseSource *> selectedSources;
    for (const SbString &requested : transaction.sourceInstanceKeys) {
	const std::string key = normalized_index_key(requested.getString());
	if (key.empty())
	    return -1;
	if (!requestedSourceKeys.insert(key).second)
	    continue;
	SoBRLDatabaseSource *source =
	    this->findIndexedDatabaseSourceInstance(key.c_str());
	if (!source || !selectedSources.insert(source).second)
	    continue;
	SoGroup *parent =
	    this->findIndexedDatabaseSourceInstanceParent(key.c_str());
	if (!parent || scene_group_find_child_index(parent, source) < 0)
	    continue;
	BObolDatabaseSourceSummary summary;
	if (!this->databaseSourceSummaryForSource(source, summary))
	    continue;
	const std::string parentPath = normalized_index_key(
	    summary.parentGroupPath.getString());
	const bool covered = std::any_of(roots.begin(), roots.end(),
	    [&parentPath](const GroupTarget &selected) {
		return scene_path_contains(selected.path, parentPath);
	    });
	if (!covered)
	    sources.push_back({source, parent});
    }

    SceneChildOrders orders;
    const auto removeEdge = [&orders](SoGroup *parent, SoNode *target) {
	auto inserted = orders.emplace(parent, std::vector<SoNode *>());
	auto &children = inserted.first->second;
	if (inserted.second) {
	    children.reserve(size_t(parent->getNumChildren()));
	    for (int i = 0; i < parent->getNumChildren(); ++i)
		children.push_back(parent->getChild(i));
	}
	auto found = std::find(children.begin(), children.end(), target);
	if (found == children.end())
	    return false;
	children.erase(found);
	return true;
    };

    int removed = 0;
    for (const GroupTarget &target : roots)
	removed += removeEdge(target.parent, target.group) ? 1 : 0;
    for (const SourceTarget &target : sources)
	removed += removeEdge(target.parent, target.source) ? 1 : 0;
    if (!removed)
	return 0;

    ChildPublication publication(*this, orders, true);
    return publication.publish();
}

int
BObolSceneController::clearDatabaseSources(void)
{
    if (!this->d->root || !this->d->root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(this->d->root);
    if (!group->getNumChildren()) return 0;
    return ChildPublication::clearSources(*this, *group);
}

static int
scene_path_component_count(const char *path);

static SbString
database_source_parent_path_from_source_path(const char *sourcePath)
{
    const char *start = skip_leading_slash(sourcePath);
    if (!start || !start[0])
	return "/";

    size_t len = strlen(start);
    while (len > 0 && start[len - 1] == '/')
	len--;
    if (len == 0)
	return "/";

    const char *slash = NULL;
    for (size_t i = 0; i < len; i++) {
	if (start[i] == '/')
	    slash = &start[i];
    }
    if (!slash)
	return "/";

    return SbString(std::string(start, (size_t)(slash - start)).c_str());
}

SbBool
BObolSceneController::databaseSourceSummaryForSource(
    SoBRLDatabaseSource *source,
    BObolDatabaseSourceSummary &summary) const
{
    summary = BObolDatabaseSourceSummary();
    if (!source || !source->getSummary(summary) || !summary.valid)
	return FALSE;

    const SbString instanceKey = database_source_effective_instance_key(source);
    SoGroup *parent = this->findIndexedDatabaseSourceInstanceParent(
			  instanceKey.getString());
    SbString parentGroupPath("");
    if (parent) {
	parentGroupPath = scene_group_summary_path(parent, "");
	if (parentGroupPath.getLength() == 0 &&
	    parent == scene_root_group(this->d->root))
	    parentGroupPath = "/";
    }
    if (parentGroupPath.getLength() == 0)
	parentGroupPath = database_source_parent_path_from_source_path(
			      source->path.getValue().getString());

    summary.hasParent = TRUE;
    summary.parentGroupPath = parentGroupPath;
    summary.drawTreeDepth =
	scene_path_component_count(parentGroupPath.getString()) + 1;
    return TRUE;
}

SbBool
BObolSceneController::getDatabaseSourceSummary(int index,
	BObolDatabaseSourceSummary &summary) const
{
    summary = BObolDatabaseSourceSummary();
    if (index < 0 || !this->d->root ||
	!this->d->root->isOfType(SoGroup::getClassTypeId()))
	return FALSE;

    return this->databaseSourceSummaryForSource(
	       this->getDatabaseSource(index), summary);
}

SbBool
BObolSceneController::getDatabaseSourceSummaryForPath(
    const char *sourcePath,
    BObolDatabaseSourceSummary &summary) const
{
    summary = BObolDatabaseSourceSummary();
    if (!sourcePath || !sourcePath[0])
	return FALSE;
    return this->databaseSourceSummaryForSource(
	       this->findIndexedDatabaseSource(sourcePath), summary);
}

int
BObolSceneController::getDatabaseSourceInstanceCountForPath(
    const char *sourcePath) const
{
    if (!sourcePath || !sourcePath[0])
	return 0;
    if (!this->d->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();

    auto it = this->d->indexes.databaseSourcePathInstancesIndex.find(
	normalized_index_key(sourcePath));
    return it == this->d->indexes.databaseSourcePathInstancesIndex.end() ? 0 :
	static_cast<int>(it->second.size());
}

SbBool
BObolSceneController::getDatabaseSourceInstanceSummaryForPath(
    const char *sourcePath,
    int instanceIndex,
    BObolDatabaseSourceSummary &summary) const
{
    summary = BObolDatabaseSourceSummary();
    if (!sourcePath || !sourcePath[0] || instanceIndex < 0)
	return FALSE;
    if (!this->d->databaseSourceIndexValid)
	this->rebuildDatabaseSourceIndex();

    auto it = this->d->indexes.databaseSourcePathInstancesIndex.find(
	normalized_index_key(sourcePath));
    if (it == this->d->indexes.databaseSourcePathInstancesIndex.end() ||
	static_cast<size_t>(instanceIndex) >= it->second.size())
	return FALSE;
    return this->databaseSourceSummaryForSource(
	it->second[static_cast<size_t>(instanceIndex)], summary);
}

SbBool
BObolSceneController::getDatabaseSourceSummaryForInstance(
    const char *sourceInstanceKey,
    BObolDatabaseSourceSummary &summary) const
{
    summary = BObolDatabaseSourceSummary();
    if (!sourceInstanceKey || !sourceInstanceKey[0])
	return FALSE;
    return this->databaseSourceSummaryForSource(
	       this->findIndexedDatabaseSourceInstance(sourceInstanceKey),
	       summary);
}

int
BObolSceneController::getRealizedShapeSummaryCount(void) const
{
    int count = 0;
    for (int i = 0; i < this->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = this->getDatabaseSource(i);
	if (source)
	    count += source->getRealizedShapeSummaryCount();
    }
    return count;
}

SbBool
BObolSceneController::getRealizedShapeSummary(int index,
	BObolRealizedShapeSummary &summary) const
{
    summary = BObolRealizedShapeSummary();
    if (index < 0)
	return FALSE;

    int remaining = index;
    for (int i = 0; i < this->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = this->getDatabaseSource(i);
	if (!source)
	    continue;

	const int sourceShapeCount =
	    source->getRealizedShapeSummaryCount();
	if (remaining < sourceShapeCount) {
	    if (!source->getRealizedShapeSummary(remaining, summary))
		return FALSE;
	    summary.ownerSourceIndex = i;
	    if (summary.ownerSourceInstanceKey.getLength() == 0)
		summary.ownerSourceInstanceKey =
		    database_source_effective_instance_key(source);
	    return TRUE;
	}
	remaining -= sourceShapeCount;
    }

    return FALSE;
}

int
BObolSceneController::getRealizedMaterialSummaryCount(void) const
{
    int count = 0;
    for (int i = 0; i < this->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = this->getDatabaseSource(i);
	if (source)
	    count += source->getRealizedMaterialSummaryCount();
    }
    return count;
}

SbBool
BObolSceneController::getRealizedMaterialSummary(int index,
	BObolRealizedMaterialSummary &summary) const
{
    summary = BObolRealizedMaterialSummary();
    if (index < 0)
	return FALSE;

    int remaining = index;
    for (int i = 0; i < this->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = this->getDatabaseSource(i);
	if (!source)
	    continue;

	const int sourceMaterialCount =
	    source->getRealizedMaterialSummaryCount();
	if (remaining < sourceMaterialCount) {
	    if (!source->getRealizedMaterialSummary(remaining, summary))
		return FALSE;
	    summary.ownerSourceIndex = i;
	    if (summary.ownerSourceInstanceKey.getLength() == 0)
		summary.ownerSourceInstanceKey =
		    database_source_effective_instance_key(source);
	    return TRUE;
	}
	remaining -= sourceMaterialCount;
    }

    return FALSE;
}

SbBool
BObolSceneController::getRealizedMaterialProperty(int materialIndex,
	int propertyIndex, SbString &groupOut, SbString &nameOut,
	SbString &valueOut) const
{
    if (materialIndex < 0 || propertyIndex < 0)
	return FALSE;

    int remaining = materialIndex;
    for (int i = 0; i < this->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = this->getDatabaseSource(i);
	if (!source)
	    continue;

	const int sourceMaterialCount =
	    source->getRealizedMaterialSummaryCount();
	if (remaining < sourceMaterialCount)
	    return source->getRealizedMaterialProperty(remaining,
		    propertyIndex, groupOut, nameOut, valueOut);
	remaining -= sourceMaterialCount;
    }

    return FALSE;
}

static void
scene_tree_summary_fill(const SoNode *node, int depth, SbBool hasParent,
			int ownerSourceIndex, const SbString &ownerSourcePath,
			const SbString &nodePath, BObolSceneTreeSummary &summary)
{
    summary = BObolSceneTreeSummary();
    if (!node)
	return;

    summary.valid = TRUE;
    summary.hasParent = hasParent;
    summary.drawTreeDepth = depth;
    summary.ownerSourceIndex = ownerSourceIndex;
    summary.ownerSourcePath = ownerSourcePath;
    if (summary.ownerSourcePath.getLength() == 0) {
	const SbString retainedOwnerPath = scene_shape_owner_source_path(node);
	if (retainedOwnerPath.getLength() > 0)
	    summary.ownerSourcePath = retainedOwnerPath;
    }
    summary.ownerSourceInstanceKey =
	scene_shape_owner_source_instance_key(node);
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
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	summary.nodeKind = BObolSceneTreeSummary::NODE_DATABASE_SOURCE;
	summary.path = source->path.getValue();
	summary.displayName = source->displayName.getValue();
	if (summary.ownerSourcePath.getLength() == 0)
	    summary.ownerSourcePath = source->path.getValue();
	if (summary.ownerSourceInstanceKey.getLength() == 0)
	    summary.ownerSourceInstanceKey =
		database_source_effective_instance_key(source);
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
    summary.path = scene_group_summary_path(node, nodePath);
}


static const SoNode *
scene_public_realized_shape_node(const SoNode *node)
{
    if (!node)
	return NULL;

    if (node->isOfType(SoBRLVListShape::getClassTypeId()) ||
	node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return node;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    const SoNode *found =
		scene_public_realized_shape_node(group->getChild(i));
	    if (found)
		return found;
	}
    }

    return NULL;
}


static const SoNode *
scene_public_realized_material_node(const SoNode *node)
{
    if (!node)
	return NULL;

    if (node->isOfType(SoBRLMaterialObject::getClassTypeId()))
	return node;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    const SoNode *found =
		scene_public_realized_material_node(group->getChild(i));
	    if (found)
		return found;
	}
    }

    return NULL;
}


static const SoNode *
scene_public_realized_child_node(const SoNode *node)
{
    const SoNode *shape = scene_public_realized_shape_node(node);
    if (shape)
	return shape;

    const SoNode *material = scene_public_realized_material_node(node);
    if (material)
	return material;

    return node;
}


static int
scene_tree_summary_node_count(const SoNode *node)
{
    if (!node)
	return 0;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	return source->getRealizedTreeSummaryCount();
    }

    int count = 1;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    count += scene_tree_summary_node_count(group->getChild(i));
    }
    return count;
}

static SbBool
find_scene_tree_summary_in_node(const SoNode *node, int &index, int depth,
				SbBool hasParent, int ownerSourceIndex,
				const SbString &ownerSourcePath, const SbString &nodePath,
				BObolSceneTreeSummary &summary)
{
    if (!node)
	return FALSE;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	const int sourceTreeCount = source->getRealizedTreeSummaryCount();
	if (index >= sourceTreeCount) {
	    index -= sourceTreeCount;
	    return FALSE;
	}
	const int sourceTreeIndex = index;
	if (!source->getRealizedTreeSummary(sourceTreeIndex, summary))
	    return FALSE;
	summary.drawTreeDepth += depth;
	summary.hasParent = hasParent ? TRUE : summary.hasParent;
	summary.ownerSourceIndex = ownerSourceIndex;
	if (ownerSourcePath.getLength() > 0)
	    summary.ownerSourcePath = ownerSourcePath;
	if (summary.ownerSourceInstanceKey.getLength() == 0)
	    summary.ownerSourceInstanceKey =
		database_source_effective_instance_key(source);
	return TRUE;
    }

    if (index == 0) {
	scene_tree_summary_fill(node, depth, hasParent, ownerSourceIndex,
				ownerSourcePath, nodePath, summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    const SoNode *child = group->getChild(i);
	    const SbString childPath =
		scene_child_summary_path(nodePath, child);
	    if (find_scene_tree_summary_in_node(child, index,
						depth + 1, TRUE, ownerSourceIndex, ownerSourcePath,
						childPath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
BObolSceneController::getSceneTreeSummaryCount(void) const
{
    return scene_tree_summary_node_count(this->d->root);
}

SbBool
BObolSceneController::getSceneTreeSummary(int index,
	BObolSceneTreeSummary &summary) const
{
    summary = BObolSceneTreeSummary();
    if (index < 0 || !this->d->root)
	return FALSE;

    if (index == 0) {
	scene_tree_summary_fill(this->d->root, 0, FALSE, -1, "", "/",
				summary);
	return TRUE;
    }
    index--;

    if (!this->d->root->isOfType(SoGroup::getClassTypeId()))
	return FALSE;

    const SoGroup *group = static_cast<const SoGroup *>(this->d->root);
    int sourceIndex = 0;
    for (int i = 0; i < group->getNumChildren(); i++) {
	const SoNode *child = group->getChild(i);
	int childOwnerIndex = -1;
	SbString childOwnerPath("");
	if (child &&
	    child->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    const SoBRLDatabaseSource *source =
		static_cast<const SoBRLDatabaseSource *>(child);
	    childOwnerIndex = sourceIndex;
	    childOwnerPath = source->path.getValue();
	    sourceIndex++;
	}

	const SbString childPath = scene_child_summary_path("/", child);
	if (find_scene_tree_summary_in_node(child, index, 1, TRUE,
					    childOwnerIndex, childOwnerPath, childPath, summary))
	    return TRUE;
    }

    return FALSE;
}

static int
scene_summary_path_equal(const SbString &summaryPath, const char *nodePath)
{
    return scene_path_equal(summaryPath.getString(), nodePath);
}

static int
scene_path_component_count(const char *path)
{
    if (!path || !path[0])
	return 0;

    int count = 0;
    const char *cp = path;
    while (*cp) {
	while (*cp == '/')
	    cp++;
	if (!*cp)
	    break;
	count++;
	while (*cp && *cp != '/')
	    cp++;
    }

    return count;
}

static const SoGroup *
scene_tree_group_find_by_summary_path(const SoNode *node,
				      const char *nodePath,
				      const SbString &fallbackPath)
{
    if (!node || !nodePath || !node->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    const SoGroup *group = static_cast<const SoGroup *>(node);
    const SbString groupPath = scene_group_summary_path(node, fallbackPath);
    if (scene_path_equal(groupPath.getString(), nodePath))
	return group;

    for (int i = 0; i < group->getNumChildren(); i++) {
	const SoNode *child = group->getChild(i);
	if (!child || !child->isOfType(SoGroup::getClassTypeId()))
	    continue;
	const SbString childPath = scene_child_summary_path(groupPath, child);
	const SoGroup *found = scene_tree_group_find_by_summary_path(child,
			       nodePath, childPath);
	if (found)
	    return found;
    }

    return NULL;
}

SbBool
BObolSceneController::getSceneTreeSummaryForPath(const char *nodePath,
	BObolSceneTreeSummary &summary) const
{
    summary = BObolSceneTreeSummary();
    if (!this->d->root || !nodePath)
	return FALSE;

    const char *normalizedPath = skip_leading_slash(nodePath);
    if (!nodePath[0] || !normalizedPath[0]) {
	scene_tree_summary_fill(this->d->root, 0, FALSE, -1, "", "/",
				summary);
	return summary.valid;
    }

    SoBRLDatabaseSource *source = this->findIndexedDatabaseSource(nodePath);
    if (source) {
	const int depth = scene_path_component_count(nodePath);
	scene_tree_summary_fill(source, depth, depth > 0 ? TRUE : FALSE, -1,
				source->path.getValue(), source->path.getValue(), summary);
	return summary.valid;
    }

    const SoGroup *group = scene_group_find_path_const(this->d->root, nodePath);
    if (group) {
	const int depth = scene_path_component_count(nodePath);
	scene_tree_summary_fill(group, depth, depth > 0 ? TRUE : FALSE, -1,
				"", nodePath, summary);
	return summary.valid;
    }

    const SoNode *shape = scene_shape_find_path(this->d->root, nodePath, NULL);
    if (shape) {
	const int depth = scene_path_component_count(nodePath);
	scene_tree_summary_fill(shape, depth, depth > 0 ? TRUE : FALSE, -1,
				"", nodePath, summary);
	return summary.valid;
    }

    if (this->getCompactSceneTreeSummaryForPath(nodePath, FALSE, summary))
	return TRUE;

    BObolSceneTreeSummary candidate;
    BObolSceneTreeSummary fallback;
    const int count = this->getSceneTreeSummaryCount();
    for (int i = 0; i < count; i++) {
	if (!this->getSceneTreeSummary(i, candidate) ||
	    !candidate.valid ||
	    !scene_summary_path_equal(candidate.path, nodePath))
	    continue;
	if (candidate.nodeKind ==
	    BObolSceneTreeSummary::NODE_DATABASE_SOURCE) {
	    summary = candidate;
	    return TRUE;
	}
	if (!fallback.valid)
	    fallback = candidate;
    }

    if (fallback.valid) {
	summary = fallback;
	return TRUE;
    }

    return FALSE;
}

SbBool
BObolSceneController::getCompactSceneTreeSummaryForPath(
    const char *nodePath, SbBool includeDescendants,
    BObolSceneTreeSummary &summary) const
{
    summary = BObolSceneTreeSummary();
    if (!this->d->root || !nodePath || !nodePath[0])
	return FALSE;

    /* Compact occurrences deliberately have no per-leaf Coin node, but they
     * are still semantic scene entities that selection and edit clients must
     * be able to address.  Resolve one through each source's ordered compact
     * path index without materializing geometry or a permanent scene node. */
    for (int sourceIndex = 0; sourceIndex < this->getDatabaseSourceCount();
	sourceIndex++) {
	SoBRLDatabaseSource *compactSource =
	    this->getDatabaseSource(sourceIndex);
	BObolCompactInstanceHandle compactHandle;
	BObolCompactInstanceSummary compactSummary;
	if (!compactSource || !compactSource->visible.getValue() ||
	    !compactSource->getCompactInstanceForPath(nodePath,
		includeDescendants, TRUE, compactHandle, compactSummary) ||
	    !compactSummary.valid || !compactSummary.visible)
	    continue;

	summary.valid = TRUE;
	summary.nodeKind = BObolSceneTreeSummary::NODE_DATABASE_SOURCE;
	summary.isDatabaseSource = TRUE;
	summary.hasParent = TRUE;
	summary.drawTreeDepth = scene_path_component_count(
	    compactSummary.path.getString());
	summary.ownerSourceIndex = sourceIndex;
	summary.ownerSourcePath = compactSource->path.getValue();
	summary.ownerSourceInstanceKey =
	    database_source_effective_instance_key(compactSource);
	summary.path = compactSummary.path;
	summary.sourceName = compactSummary.sourceName;
	summary.sourceType = compactSummary.geometryKind;
	summary.displayName = compactSummary.sourceName;
	summary.geometryName = compactSummary.sourceName;
	return TRUE;
    }

    return FALSE;
}

SbBool
BObolSceneController::getSceneChildTreeSummary(const char *nodePath,
	int childIndex,
	BObolSceneTreeSummary &summary) const
{
    summary = BObolSceneTreeSummary();
    if (!this->d->root || !nodePath || childIndex < 0)
	return FALSE;

    BObolSceneTreeSummary parentSummary;
    if (!this->getSceneTreeSummaryForPath(nodePath, parentSummary))
	return FALSE;

    SoBRLDatabaseSource *source = this->findDatabaseSource(nodePath);
    if (source) {
	if (childIndex >= source->getNumChildren())
	    return FALSE;

	const SoNode *child = source->getChild(childIndex);
	const SoNode *publicChild =
	    (child && child->isOfType(SoBRLDatabaseSource::getClassTypeId())) ?
	    child : scene_public_realized_child_node(child);
	const SbString childPath =
	    scene_child_summary_path(parentSummary.path, child);
	scene_tree_summary_fill(publicChild,
				parentSummary.drawTreeDepth + 1, TRUE,
				parentSummary.ownerSourceIndex,
				source->path.getValue(), childPath, summary);
	if (!summary.valid)
	    return FALSE;
	if (summary.ownerSourcePath.getLength() == 0)
	    summary.ownerSourcePath = source->path.getValue();
	if (summary.ownerSourceInstanceKey.getLength() == 0)
	    summary.ownerSourceInstanceKey =
		database_source_effective_instance_key(source);
	if (summary.path.getLength() == 0)
	    summary.path = source->path.getValue();
	return TRUE;
    }

    const SoGroup *group = this->findGroup(nodePath);
    if (!group)
	group = scene_tree_group_find_by_summary_path(this->d->root, nodePath,
		SbString("/"));
    if (group) {
	if (childIndex >= group->getNumChildren())
	    return FALSE;

	const SoNode *child = group->getChild(childIndex);
	const SbString childPath = scene_child_summary_path(
				       parentSummary.path, child);
	scene_tree_summary_fill(child, parentSummary.drawTreeDepth + 1,
				TRUE, -1, "", childPath, summary);
	return summary.valid;
    }

    return FALSE;
}

static void
scene_display_summary_fill_common(BObolSceneDisplaySummary &summary,
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
scene_display_summary_fill_shape(const ShapeT *shape, int nodeKind,
				 int ownerSourceIndex, const SbString &ownerSourcePath,
				 BObolSceneDisplaySummary &summary)
{
    SbString effectiveOwnerPath = ownerSourcePath;
    if (shape && effectiveOwnerPath.getLength() == 0 &&
	shape->ownerSourcePath.getValue().getLength() > 0)
	effectiveOwnerPath = shape->ownerSourcePath.getValue();
    scene_display_summary_fill_common(summary, nodeKind, ownerSourceIndex,
				      effectiveOwnerPath,
				      shape ? shape->sourcePath.getValue() : SbString(""));
    if (!shape)
	return;

    summary.ownerSourceInstanceKey = shape->ownerSourceInstanceKey.getValue();
    summary.hasDrawIntent = shape->sourcePath.getValue().getLength() > 0;
    summary.intentPath = shape->sourcePath.getValue();
    summary.intentDrawMode = shape->drawMode.getValue();
    summary.visible = shape->visible.getValue();
    summary.selected = shape->selected.getValue();
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
scene_display_summary_fill(const SoNode *node, int ownerSourceIndex,
			   const SbString &ownerSourcePath, const SbString &nodePath,
			   BObolSceneDisplaySummary &summary)
{
    summary = BObolSceneDisplaySummary();
    if (!node)
	return;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	scene_display_summary_fill_common(summary,
					  BObolSceneTreeSummary::NODE_DATABASE_SOURCE,
					  ownerSourceIndex, ownerSourcePath, source->path.getValue());
	summary.ownerSourceInstanceKey =
	    database_source_effective_instance_key(source);
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
	scene_display_summary_fill_shape(
	    static_cast<const SoBRLVListShape *>(node),
	    BObolSceneTreeSummary::NODE_VLIST_SHAPE,
	    ownerSourceIndex, ownerSourcePath, summary);
	return;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	scene_display_summary_fill_shape(
	    static_cast<const SoBRLMeshShape *>(node),
	    BObolSceneTreeSummary::NODE_MESH_SHAPE,
	    ownerSourceIndex, ownerSourcePath, summary);
	return;
    }

    if (node->isOfType(SoBRLMaterialObject::getClassTypeId())) {
	const SoBRLMaterialObject *object =
	    static_cast<const SoBRLMaterialObject *>(node);
	scene_display_summary_fill_common(summary,
					  BObolSceneTreeSummary::NODE_MATERIAL_OBJECT,
					  ownerSourceIndex, ownerSourcePath,
					  object->sourcePath.getValue());
	return;
    }

    if (node->isOfType(SoBRLSceneGroup::getClassTypeId())) {
	const SoBRLSceneGroup *group =
	    static_cast<const SoBRLSceneGroup *>(node);
	const SbString retainedPath =
	    scene_group_summary_path(node, nodePath);
	scene_display_summary_fill_common(summary,
					  BObolSceneTreeSummary::NODE_GROUP,
					  ownerSourceIndex, ownerSourcePath, retainedPath);
	summary.hasDrawIntent = group->drawIntentValid.getValue();
	if (summary.hasDrawIntent) {
	    summary.intentPath = group->drawIntentPath.getValue();
	    if (summary.intentPath.getLength() == 0)
		summary.intentPath = retainedPath;
	    summary.intentDrawMode = group->drawMode.getValue();
	}
	summary.visible = group->visible.getValue();
	summary.selected = group->selected.getValue();
	summary.highlighted = group->highlighted.getValue();
	summary.lineStyle = group->lineStyle.getValue();
	summary.lineWidth = group->lineWidth.getValue();
	summary.transparency = group->transparency.getValue();
	summary.drawMode = group->drawMode.getValue();
	summary.materialValid =
	    group->materialColorValid.getValue() ||
	    group->colorOverride.getValue();
	summary.materialRevision = group->materialRevision.getValue();
	if (group->materialColorValid.getValue())
	    summary.materialColor = group->materialColor.getValue();
	else if (group->colorOverride.getValue())
	    summary.materialColor = group->color.getValue();
	return;
    }

    const int nodeKind = node->isOfType(SoGroup::getClassTypeId()) ?
			 BObolSceneTreeSummary::NODE_GROUP :
			 BObolSceneTreeSummary::NODE_OTHER;
    scene_display_summary_fill_common(summary, nodeKind, ownerSourceIndex,
				      ownerSourcePath, nodePath);
}

static SbBool
find_scene_display_summary_in_node(const SoNode *node, int &index,
				   int ownerSourceIndex, const SbString &ownerSourcePath,
				   const SbString &nodePath, BObolSceneDisplaySummary &summary)
{
    if (!node)
	return FALSE;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	const int sourceDisplayCount =
	    source->getRealizedDisplaySummaryCount();
	if (index >= sourceDisplayCount) {
	    index -= sourceDisplayCount;
	    return FALSE;
	}
	const int sourceDisplayIndex = index;
	if (!source->getRealizedDisplaySummary(sourceDisplayIndex, summary))
	    return FALSE;
	summary.ownerSourceIndex = ownerSourceIndex;
	if (ownerSourcePath.getLength() > 0)
	    summary.ownerSourcePath = ownerSourcePath;
	if (summary.ownerSourceInstanceKey.getLength() == 0)
	    summary.ownerSourceInstanceKey =
		database_source_effective_instance_key(source);
	return TRUE;
    }

    if (index == 0) {
	scene_display_summary_fill(node, ownerSourceIndex, ownerSourcePath,
				   nodePath, summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    const SoNode *child = group->getChild(i);
	    const SbString childPath =
		scene_child_summary_path(nodePath, child);
	    if (find_scene_display_summary_in_node(child, index,
						   ownerSourceIndex, ownerSourcePath, childPath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
BObolSceneController::getSceneDisplaySummaryCount(void) const
{
    return this->getSceneTreeSummaryCount();
}

SbBool
BObolSceneController::getSceneDisplaySummary(int index,
	BObolSceneDisplaySummary &summary) const
{
    summary = BObolSceneDisplaySummary();
    if (index < 0 || !this->d->root)
	return FALSE;

    if (index == 0) {
	scene_display_summary_fill(this->d->root, -1, "", "/", summary);
	return TRUE;
    }
    index--;

    if (!this->d->root->isOfType(SoGroup::getClassTypeId()))
	return FALSE;

    const SoGroup *group = static_cast<const SoGroup *>(this->d->root);
    int sourceIndex = 0;
    for (int i = 0; i < group->getNumChildren(); i++) {
	const SoNode *child = group->getChild(i);
	int childOwnerIndex = -1;
	SbString childOwnerPath("");
	if (child &&
	    child->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    const SoBRLDatabaseSource *source =
		static_cast<const SoBRLDatabaseSource *>(child);
	    childOwnerIndex = sourceIndex;
	    childOwnerPath = source->path.getValue();
	    sourceIndex++;
	}

	const SbString childPath = scene_child_summary_path("/", child);
	if (find_scene_display_summary_in_node(child, index,
					       childOwnerIndex, childOwnerPath, childPath, summary))
	    return TRUE;
    }

    return FALSE;
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
BObolSceneController::getSceneMaterialSummaryCount(void) const
{
    return this->getSceneDisplaySummaryCount();
}

SbBool
BObolSceneController::getSceneMaterialSummary(int index,
	BObolSceneMaterialSummary &summary) const
{
    summary = BObolSceneMaterialSummary();
    BObolSceneDisplaySummary display;
    if (!this->getSceneDisplaySummary(index, display))
	return FALSE;

    scene_material_summary_from_display(display, summary);
    return TRUE;
}

static int
scene_bounds_node_kind(const SoNode *node)
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
scene_bounds_node_path(const SoNode *node)
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
scene_bounds_for_vlist_shape(const SoBRLVListShape *shape, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!shape)
	return FALSE;
    for (int i = 0; i < shape->point.getNum(); i++)
	bounds.extendBy(shape->point[i]);
    return shape->point.getNum() > 0;
}

static SbBool
scene_bounds_for_mesh_shape(const SoBRLMeshShape *shape, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!shape)
	return FALSE;
    for (int i = 0; i < shape->point.getNum(); i++)
	bounds.extendBy(shape->point[i]);
    return shape->point.getNum() > 0;
}

static SbBox3f
scene_bounds_transform_box(const SbBox3f &bounds, const SbMatrix &matrix)
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

static SbBool
scene_node_is_overlay_intent(const SoNode *node)
{
    if (!node)
	return FALSE;
    if (node->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return static_cast<const SoBRLSceneGroup *>(node)->
	       overlayIntent.getValue();
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<const SoBRLVListShape *>(node)->
	       overlayIntent.getValue();
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<const SoBRLMeshShape *>(node)->
	       overlayIntent.getValue();
    if (node->isOfType(SoBRLGrid::getClassTypeId()))
	return static_cast<const SoBRLGrid *>(node)->
	       overlayIntent.getValue();
    return FALSE;
}

static SbBool
scene_database_source_uses_realized_placement(
    const SoBRLDatabaseSource *source)
{
    if (!source ||
	source->realizationStatus.getValue() !=
	SoBRLDatabaseSource::REALIZED ||
	(source->realizationRoleFlags.getValue() &
	 SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL))
	return FALSE;


    return (source->hasRealizedWireGeometry() ||
	    source->hasRealizedMeshGeometry() ||
	    source->getRealizedMaterialObjectCount() > 0) ? TRUE : FALSE;
}

static SbBool
scene_bounds_for_node_transformed(const SoNode *node, const SbMatrix &matrix,
				  SbBox3f &bounds, SbBool includeOverlays)
{
    bounds.makeEmpty();
    if (!node)
	return FALSE;

    if (!includeOverlays && scene_node_is_overlay_intent(node))
	return FALSE;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SbBox3f localBounds;
	if (!scene_bounds_for_vlist_shape(
		static_cast<const SoBRLVListShape *>(node), localBounds))
	    return FALSE;
	bounds = scene_bounds_transform_box(localBounds, matrix);
	return bounds.isEmpty() ? FALSE : TRUE;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SbBox3f localBounds;
	if (!scene_bounds_for_mesh_shape(
		static_cast<const SoBRLMeshShape *>(node), localBounds))
	    return FALSE;
	bounds = scene_bounds_transform_box(localBounds, matrix);
	return bounds.isEmpty() ? FALSE : TRUE;
    }

    SbBool valid = FALSE;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	SbMatrix childMatrix = matrix;
	if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    const SoBRLDatabaseSource *source =
		static_cast<const SoBRLDatabaseSource *>(node);
	    SbBool hasSourceTransform = FALSE;
	    for (int i = 0; i < group->getNumChildren(); i++) {
		const SoNode *child = group->getChild(i);
		if (child && child->isOfType(
			SoMatrixTransform::getClassTypeId())) {
		    hasSourceTransform = TRUE;
		    break;
		}
	    }
	    if (!hasSourceTransform &&
		!scene_database_source_uses_realized_placement(source) &&
		source->drawMatrixValid.getValue())
		childMatrix.multRight(source->drawMatrix.getValue());
	}
	for (int i = 0; i < group->getNumChildren(); i++) {
	    const SoNode *child = group->getChild(i);
	    if (child &&
		child->isOfType(SoMatrixTransform::getClassTypeId())) {
		const SoMatrixTransform *transform =
		    static_cast<const SoMatrixTransform *>(child);
		childMatrix.multRight(transform->matrix.getValue());
		continue;
	    }
	    SbBox3f childBounds;
	    if (scene_bounds_for_node_transformed(child, childMatrix,
						  childBounds, includeOverlays)) {
		bounds.extendBy(childBounds);
		valid = TRUE;
	    }
	}
    }

    if (valid)
	return TRUE;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	SbBox3f localBounds;
	if (!source->getSourceBounds(localBounds))
	    return FALSE;

	SbMatrix sourceMatrix = matrix;
	if (source->drawMatrixValid.getValue())
	    sourceMatrix.multRight(source->drawMatrix.getValue());
	bounds = scene_bounds_transform_box(localBounds, sourceMatrix);
	return bounds.isEmpty() ? FALSE : TRUE;
    }

    return valid;
}

static SbBool
scene_bounds_for_node(const SoNode *node, SbBox3f &bounds,
		      SbBool includeOverlays)
{
    return scene_bounds_for_node_transformed(node, SbMatrix::identity(),
	    bounds, includeOverlays);
}

static SbBox3f
scene_autoview_padded_bounds(const SbBox3f &bounds)
{
    SbBox3f padded;
    padded.makeEmpty();
    if (bounds.isEmpty())
	return padded;

    const SbVec3f bmin = bounds.getMin();
    const SbVec3f bmax = bounds.getMax();
    const SbVec3f center = bounds.getCenter();
    float size = bmax[0] - bmin[0];
    if (bmax[1] - bmin[1] > size)
	size = bmax[1] - bmin[1];
    if (bmax[2] - bmin[2] > size)
	size = bmax[2] - bmin[2];

    const float halfSize = 0.5f * size;
    padded.extendBy(SbVec3f(center[0] - halfSize, center[1] - halfSize,
			    center[2] - halfSize));
    padded.extendBy(SbVec3f(center[0] + halfSize, center[1] + halfSize,
			    center[2] + halfSize));
    return padded;
}

static SbBool
scene_database_source_bounds_recursive(const SoNode *node,
				       SbBox3f &bounds,
				       SbBool padForAutoview)
{
    if (!node)
	return FALSE;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	SbBox3f sourceBounds;
	if (!scene_bounds_for_node(source, sourceBounds, TRUE))
	    (void)source->getEffectiveSourceBounds(sourceBounds);
	if (sourceBounds.isEmpty())
	    return FALSE;
	if (padForAutoview)
	    sourceBounds = scene_autoview_padded_bounds(sourceBounds);
	bounds.extendBy(sourceBounds);
	return sourceBounds.isEmpty() ? FALSE : TRUE;
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return FALSE;

    SbBool valid = FALSE;
    const SoGroup *group = static_cast<const SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++) {
	if (scene_database_source_bounds_recursive(group->getChild(i),
	    bounds, padForAutoview))
	    valid = TRUE;
    }

    return valid;
}

SbBool
BObolSceneController::getDatabaseSourceBounds(SbBox3f &bounds,
	SbBool padForAutoview) const
{
    bounds.makeEmpty();
    return scene_database_source_bounds_recursive(this->d->root, bounds,
	    padForAutoview);
}

static void
scene_bounds_summary_fill(const SoNode *node, int ownerSourceIndex,
			  const SbString &ownerSourcePath, const SbString &nodePath,
			  BObolSceneBoundsSummary &summary)
{
    summary = BObolSceneBoundsSummary();
    if (!node)
	return;

    summary.valid = TRUE;
    summary.nodeKind = scene_bounds_node_kind(node);
    summary.ownerSourceIndex = ownerSourceIndex;
    summary.ownerSourcePath = ownerSourcePath;
    if (summary.ownerSourcePath.getLength() == 0) {
	const SbString retainedOwnerPath = scene_shape_owner_source_path(node);
	if (retainedOwnerPath.getLength() > 0)
	    summary.ownerSourcePath = retainedOwnerPath;
    }
    summary.ownerSourceInstanceKey =
	scene_shape_owner_source_instance_key(node);
    if (summary.ownerSourceInstanceKey.getLength() == 0 &&
	node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	summary.ownerSourceInstanceKey =
	    database_source_effective_instance_key(source);
    }
    summary.path = summary.nodeKind == BObolSceneTreeSummary::NODE_GROUP ?
		   scene_group_summary_path(node, nodePath) : scene_bounds_node_path(node);
    summary.boundsValid = scene_bounds_for_node(node, summary.bounds, TRUE);
}

SbBool
BObolSceneController::getSceneSubtreeBounds(const char *nodePath,
	SbBool includeOverlays,
	SbBox3f &bounds) const
{
    bounds.makeEmpty();
    if (!this->d->root)
	return FALSE;

    const char *path = nodePath ? nodePath : "/";
    const char *normalizedPath = skip_leading_slash(path);
    const SoNode *node = NULL;

    SbBool compactBoundsValid = FALSE;
    for (int i = 0; i < this->getDatabaseSourceCount(); i++) {
	const SoBRLDatabaseSource *source = this->getDatabaseSource(i);
	if (!source || !source->hasCompactInstanceIndex() ||
	    (!includeOverlays && source->auxiliarySource.getValue()))
	    continue;
	SbBox3f sourceBounds;
	if (source->getCompactInstanceBoundsForPath(path, TRUE, sourceBounds)) {
	    bounds.extendBy(sourceBounds);
	    compactBoundsValid = TRUE;
	}
    }
    if (compactBoundsValid)
	return TRUE;

    if (!path[0] || !normalizedPath[0])
	node = this->d->root;
    if (!node)
	node = this->findDatabaseSource(path);
    if (!node)
	node = this->findGroup(path);
    if (!node)
	node = this->findShape(path);

    return scene_bounds_for_node(node, bounds, includeOverlays);
}

static SbBool
find_scene_bounds_summary_in_node(const SoNode *node, int &index,
				  int ownerSourceIndex, const SbString &ownerSourcePath,
				  const SbString &nodePath, BObolSceneBoundsSummary &summary)
{
    if (!node)
	return FALSE;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	const int sourceBoundsCount =
	    source->getRealizedBoundsSummaryCount();
	if (index >= sourceBoundsCount) {
	    index -= sourceBoundsCount;
	    return FALSE;
	}
	const int sourceBoundsIndex = index;
	if (!source->getRealizedBoundsSummary(sourceBoundsIndex, summary))
	    return FALSE;
	summary.ownerSourceIndex = ownerSourceIndex;
	if (ownerSourcePath.getLength() > 0)
	    summary.ownerSourcePath = ownerSourcePath;
	return TRUE;
    }

    if (index == 0) {
	scene_bounds_summary_fill(node, ownerSourceIndex, ownerSourcePath,
				  nodePath, summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    const SoNode *child = group->getChild(i);
	    const SbString childPath =
		scene_child_summary_path(nodePath, child);
	    if (find_scene_bounds_summary_in_node(child, index,
						  ownerSourceIndex, ownerSourcePath, childPath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
BObolSceneController::getSceneBoundsSummaryCount(void) const
{
    return this->getSceneTreeSummaryCount();
}

SbBool
BObolSceneController::getSceneBoundsSummary(int index,
	BObolSceneBoundsSummary &summary) const
{
    summary = BObolSceneBoundsSummary();
    if (index < 0 || !this->d->root)
	return FALSE;

    if (index == 0) {
	scene_bounds_summary_fill(this->d->root, -1, "", "/", summary);
	return TRUE;
    }
    index--;

    if (!this->d->root->isOfType(SoGroup::getClassTypeId()))
	return FALSE;

    const SoGroup *group = static_cast<const SoGroup *>(this->d->root);
    int sourceIndex = 0;
    for (int i = 0; i < group->getNumChildren(); i++) {
	const SoNode *child = group->getChild(i);
	int childOwnerIndex = -1;
	SbString childOwnerPath("");
	if (child &&
	    child->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	    const SoBRLDatabaseSource *source =
		static_cast<const SoBRLDatabaseSource *>(child);
	    childOwnerIndex = sourceIndex;
	    childOwnerPath = source->path.getValue();
	    sourceIndex++;
	}

	const SbString childPath = scene_child_summary_path("/", child);
	if (find_scene_bounds_summary_in_node(child, index,
					      childOwnerIndex, childOwnerPath, childPath, summary))
	    return TRUE;
    }

    return FALSE;
}

unsigned int
BObolSceneController::getLastVisitedSourceCount(void) const
{
    return this->d->activeRealizationAction ?
	this->d->activeRealizationAction->getVisitedSourceCount() : this->d->lastVisitedSourceCount;
}

unsigned int
BObolSceneController::getLastRealizedSourceCount(void) const
{
    return this->d->activeRealizationAction ?
	this->d->activeRealizationAction->getRealizedSourceCount() : this->d->lastRealizedSourceCount;
}

unsigned int
BObolSceneController::getLastFailedSourceCount(void) const
{
    return this->d->activeRealizationAction ?
	this->d->activeRealizationAction->getFailedSourceCount() : this->d->lastFailedSourceCount;
}

const SbString &
BObolSceneController::getLastDiagnostics(void) const
{
    return this->d->activeRealizationAction ?
	this->d->activeRealizationAction->getDiagnostics() : this->d->lastDiagnostics;
}

void
BObolSceneController::advanceFrameRevision(void)
{
    if (this->d->mutationBatchDepth > 0) {
	this->d->mutationBatchFrameRevisionPending = TRUE;
	return;
    }

    bobol_identity_advance(this->d->frameRevision);
}

void
BObolSceneController::advanceStructuralRevision(SbBool stopRealizationTraversal)
{
    if (stopRealizationTraversal && this->d->activeRealizationAction)
	this->d->activeRealizationAction->stopSceneTraversal();
    if (this->d->mutationBatchDepth > 0) {
	this->d->mutationBatchStructuralRevisionPending = TRUE;
	this->d->mutationBatchFrameRevisionPending = TRUE;
	return;
    }

    bobol_identity_advance(this->d->structuralRevision);
    this->advanceFrameRevision();
}
