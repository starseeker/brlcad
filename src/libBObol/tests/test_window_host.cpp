/*              T E S T _ W I N D O W _ H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/app.h"

#include "bu/str.h"

#include "BObol.h"
#include "BObol/BImageSource.h"
#include "BObol/BViewportImage.h"

#include "bu/log.h"
#include "bu/malloc.h"

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wshadow"
#endif
#include "imgstream/fb_compat.h"
#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

#include <Inventor/SbViewportRegion.h>
#include <Inventor/nodes/SoClipPlane.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoTexture2.h>

#include <math.h>
#include <string.h>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <future>
#include <mutex>
#include <string>
#include <thread>

#define CHECK(_expr, _msg) do { \
    if (!(_expr)) { \
	bu_log("FAIL: %s\n", _msg); \
	return 1; \
    } \
} while (0)

static int
float_equal(float a, float b)
{
    return fabs((double)a - (double)b) <= 1.0e-6;
}

static SoClipPlane *
find_clip_plane(SoNode *node, const char *name)
{
    if (!node || !name)
	return NULL;
    if (node->isOfType(SoClipPlane::getClassTypeId()) &&
	bu_strcmp(node->getName().getString(), name) == 0)
	return static_cast<SoClipPlane *>(node);
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return NULL;
    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoClipPlane *plane = find_clip_plane(group->getChild(i), name);
	if (plane)
	    return plane;
    }
    return NULL;
}

static void
property_value_init(struct bv_display_property_value *value)
{
    memset(value, 0, sizeof(*value));
    value->struct_size = sizeof(*value);
    value->type = BV_DISPLAY_PROPERTY_BOOL;
}

struct TestPropertyProvider {
    int zclip;
    double perspective;
    int framebufferMode;
    int adcVisible;
    int modelAxesVisible;
    int viewAxesVisible;
};

static int
test_property_provider_get(void *data, const char *name,
	struct bv_display_property_value *value)
{
    TestPropertyProvider *provider =
	static_cast<TestPropertyProvider *>(data);
    if (!provider || !name || !value)
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;
    if (bu_strcmp(name, "view.zclip") == 0) {
	value->bool_value = provider->zclip;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (bu_strcmp(name, "view.perspective") == 0) {
	value->double_value = provider->perspective;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (bu_strcmp(name, "composition.framebuffer.mode") == 0) {
	if (provider->framebufferMode == 0)
	    value->string_value = "off";
	else if (provider->framebufferMode == 1)
	    value->string_value = "overlay";
	else if (provider->framebufferMode == 2)
	    value->string_value = "underlay";
	else if (provider->framebufferMode == 3)
	    value->string_value = "interlay";
	else
	    return BV_DISPLAY_PROPERTY_INVALID;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (bu_strcmp(name, "view.faceplate.adc.visible") == 0) {
	value->bool_value = provider->adcVisible;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (bu_strcmp(name, "view.faceplate.model_axes.visible") == 0) {
	value->bool_value = provider->modelAxesVisible;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (bu_strcmp(name, "view.faceplate.view_axes.visible") == 0) {
	value->bool_value = provider->viewAxesVisible;
	return BV_DISPLAY_PROPERTY_OK;
    }
    return BV_DISPLAY_PROPERTY_UNSUPPORTED;
}

static int
test_property_provider_set(void *data, const char *name,
	const struct bv_display_property_value *value)
{
    TestPropertyProvider *provider =
	static_cast<TestPropertyProvider *>(data);
    if (!provider || !name || !value)
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;
    if (bu_strcmp(name, "view.zclip") == 0) {
	provider->zclip = value->bool_value ? 1 : 0;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (bu_strcmp(name, "view.perspective") == 0) {
	provider->perspective = value->double_value;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (bu_strcmp(name, "composition.framebuffer.mode") == 0) {
	if (!value->string_value)
	    return BV_DISPLAY_PROPERTY_INVALID;
	if (bu_strcmp(value->string_value, "off") == 0)
	    provider->framebufferMode = 0;
	else if (bu_strcmp(value->string_value, "overlay") == 0)
	    provider->framebufferMode = 1;
	else if (bu_strcmp(value->string_value, "underlay") == 0)
	    provider->framebufferMode = 2;
	else if (bu_strcmp(value->string_value, "interlay") == 0)
	    provider->framebufferMode = 3;
	else
	    return BV_DISPLAY_PROPERTY_INVALID;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (bu_strcmp(name, "view.faceplate.adc.visible") == 0) {
	provider->adcVisible = value->bool_value ? 1 : 0;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (bu_strcmp(name, "view.faceplate.model_axes.visible") == 0) {
	provider->modelAxesVisible = value->bool_value ? 1 : 0;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (bu_strcmp(name, "view.faceplate.view_axes.visible") == 0) {
	provider->viewAxesVisible = value->bool_value ? 1 : 0;
	return BV_DISPLAY_PROPERTY_OK;
    }
    return BV_DISPLAY_PROPERTY_UNSUPPORTED;
}

struct InputActionState {
    int calls = 0;
    BObolInputAction action = BOBOL_ACTION_NONE;
    int result = BOBOL_INPUT_RESULT_HANDLED;
};

/* Test-only action vocabulary.  These values intentionally overlap standard
 * view actions to prove dispatch is selected by binding-layer provenance. */
enum TestLayerAction {
    TEST_LAYER_SELECT_COMMIT = BOBOL_ACTION_VIEW_PRIMARY_RELEASE,
    TEST_LAYER_ROTATE_BEGIN = BOBOL_ACTION_VIEW_ROTATE,
    TEST_LAYER_ROTATE_CANCEL = BOBOL_ACTION_VIEW_PAN,
    TEST_LAYER_TRANSLATE_BEGIN = BOBOL_ACTION_VIEW_ZOOM,
    TEST_LAYER_TRANSLATE_CANCEL = BOBOL_ACTION_VIEW_CENTER,
    TEST_LAYER_SCALE_BEGIN = BOBOL_ACTION_VIEW_FRONT,
    TEST_LAYER_SCALE_CANCEL = BOBOL_ACTION_VIEW_TOP,
    TEST_LAYER_CONSTRAINED_TRANSLATE_BEGIN = BOBOL_ACTION_VIEW_BOTTOM,
    TEST_LAYER_CONSTRAINED_TRANSLATE_CANCEL = BOBOL_ACTION_VIEW_LEFT,
    TEST_LAYER_CONSTRAINED_ROTATE_BEGIN = BOBOL_ACTION_VIEW_REAR,
    TEST_LAYER_CONSTRAINED_ROTATE_CANCEL = BOBOL_ACTION_VIEW_RIGHT,
    TEST_LAYER_CONSTRAINED_SCALE_BEGIN = BOBOL_ACTION_VIEW_2,
    TEST_LAYER_CONSTRAINED_SCALE_CANCEL = BOBOL_ACTION_VIEW_3,
    TEST_LAYER_PERSPECTIVE_MODE_TOGGLE = BOBOL_ACTION_VIEW_4,
    TEST_LAYER_PERSPECTIVE_CYCLE = BOBOL_ACTION_VIEW_5,
    TEST_LAYER_FACEPLATE_TOGGLE = BOBOL_ACTION_TOGGLE_ADC,
    TEST_LAYER_FACEPLATE_GUI_TOGGLE = BOBOL_ACTION_TOGGLE_MODEL_AXES
};

static int
test_input_action(void *data, BObolInputAction action,
	const BObolInputEvent *event)
{
    InputActionState *state = static_cast<InputActionState *>(data);
    if (!state || !event)
	return -1;
    state->calls++;
    state->action = action;
    return state->result;
}

struct PendingProgressState {
    int calls = 0;
};

static int
test_pending_progress_provider(BObolViewController *controller,
	void *data, const BObolProgressiveOptions *options,
	BObolProgressiveStatus *status)
{
    PendingProgressState *state = static_cast<PendingProgressState *>(data);
    if (!controller || !state || !options || !status)
	return -1;
    state->calls++;
    status->hasMore = 1;
    return 1;
}

static int
test_input_context(void)
{
    InputActionState firstState;
    InputActionState secondState;
    BObolInputContext first;
    BObolInputContext second;
    first.setProfile(&BObolInputContext::defaultViewProfile());
    first.setActionHandler(test_input_action, &firstState);
    second.setActionHandler(test_input_action, &secondState);

    BObolInputEvent event;
    event.type = BOBOL_INPUT_POINTER_MOTION;
    event.button = 0;
    event.modifiers = BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_SHIFT;
    CHECK(first.dispatch(&event) == 1 && firstState.calls == 1 &&
	  firstState.action == BOBOL_ACTION_VIEW_ZOOM,
	  "input context chooses the most-specific modifier binding");
    CHECK(second.dispatch(&event) == 0 && secondState.calls == 0,
	  "input contexts do not share bindings or action handlers");

    event.type = BOBOL_INPUT_KEY_PRESS;
    event.key = 'R';
    event.button = BOBOL_INPUT_ANY;
    event.modifiers = BOBOL_INPUT_MOD_SHIFT;
    CHECK(first.dispatch(&event) == 1 &&
	  firstState.action == BOBOL_ACTION_VIEW_REAR,
	  "input context preserves modifier-specific view actions");
    event.type = BOBOL_INPUT_WHEEL;
    event.wheelDelta = -15;
    event.button = BOBOL_INPUT_ANY;
    event.modifiers = BOBOL_INPUT_MOD_NONE;
    CHECK(first.dispatch(&event) == 1 &&
	  firstState.action == BOBOL_ACTION_VIEW_ZOOM,
	  "input context dispatches normalized wheel zoom");
    const int actionCalls = firstState.calls;
    event.type = BOBOL_INPUT_KEY_PRESS;
    event.key = 'A';
    event.modifiers = BOBOL_INPUT_MOD_SHIFT;
    CHECK(first.dispatch(&event) == 0 && firstState.calls == actionCalls,
	  "modified application shortcuts do not trigger plain view actions");
    event.key = 'R';
    event.modifiers = BOBOL_INPUT_MOD_SHIFT | BOBOL_INPUT_MOD_CONTROL;
    CHECK(first.dispatch(&event) == 0 && firstState.calls == actionCalls,
	  "rear view requires exactly the shift modifier");
    event.type = BOBOL_INPUT_POINTER_RELEASE;
    event.button = 0;
    event.modifiers = BOBOL_INPUT_MOD_NONE;
    CHECK(first.dispatch(&event) == 1,
	  "input context dispatches pointer-release actions");
    CHECK(firstState.action == BOBOL_ACTION_VIEW_PRIMARY_RELEASE,
	  "input context normalizes pointer-release actions");
    event.type = BOBOL_INPUT_POINTER_PRESS;
    event.button = 0;
    event.modifiers = BOBOL_INPUT_MOD_SHIFT;
    CHECK(first.dispatch(&event) == 1,
	  "input context dispatches primary pan-begin actions");
    CHECK(firstState.action == BOBOL_ACTION_VIEW_PAN_BEGIN,
	  "input context normalizes primary pan-begin actions");
    event.type = BOBOL_INPUT_POINTER_RELEASE;
    CHECK(first.dispatch(&event) == 1,
	  "input context dispatches primary pan-end actions");
    CHECK(firstState.action == BOBOL_ACTION_VIEW_PAN_END,
	  "input context normalizes primary pan-end actions");
    CHECK(first.hasAction(BOBOL_ACTION_VIEW_ROTATE) &&
	  first.hasAction(BOBOL_ACTION_VIEW_PAN_BEGIN) &&
	  first.hasAction(BOBOL_ACTION_VIEW_PAN_END) &&
	  !second.hasAction(BOBOL_ACTION_VIEW_ROTATE),
	  "input contexts report their own action vocabulary");

    BObolInputBinding semantic_binding;
    semantic_binding.eventType = BOBOL_INPUT_POINTER_RELEASE;
    semantic_binding.key = BOBOL_INPUT_ANY;
    semantic_binding.button = 0;
    semantic_binding.requiredModifiers = BOBOL_INPUT_MOD_NONE;
    semantic_binding.forbiddenModifiers = BOBOL_INPUT_MOD_SHIFT |
	BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_ALT |
	BOBOL_INPUT_MOD_META;
    semantic_binding.priority = 0;
    semantic_binding.action = TEST_LAYER_SELECT_COMMIT;
    InputActionState semantic_state;
    semantic_state.result = BOBOL_INPUT_RESULT_DEFERRED;
    const BObolInputActionLayer semantic_layer = {"test-selection",
	&semantic_binding, 1, test_input_action};
    BObolInputBinding background_binding = semantic_binding;
    background_binding.priority = -1;
    const BObolInputActionLayer background_layer = {"test-background",
	&background_binding, 1, test_input_action};
    CHECK(first.setActionLayer(&semantic_layer, &semantic_state,
	  &semantic_state) &&
	  first.setActionLayer(&background_layer, &secondState, &secondState) &&
	  first.clearActionLayerIf(&secondState),
	  "application action layers coexist and tear down by guarded owner");
    event.type = BOBOL_INPUT_POINTER_RELEASE;
    event.button = 0;
    event.modifiers = BOBOL_INPUT_MOD_NONE;
    const int view_calls_before_layer = firstState.calls;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_DEFERRED &&
	  semantic_state.calls == 1 &&
	  semantic_state.action == TEST_LAYER_SELECT_COMMIT &&
	  firstState.calls == view_calls_before_layer,
	  "application layers win equal-score view bindings even when action IDs overlap");
    semantic_state.result = BOBOL_INPUT_RESULT_UNHANDLED;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_HANDLED &&
	  semantic_state.calls == 2 &&
	  firstState.calls == view_calls_before_layer + 1 &&
	  firstState.action == BOBOL_ACTION_VIEW_PRIMARY_RELEASE,
	  "unhandled application actions fall through to the matching view action");
    CHECK(!first.clearActionLayerIf(&secondState) &&
	  first.clearActionLayerIf(&semantic_state) &&
	  first.dispatch(&event) == BOBOL_INPUT_RESULT_HANDLED &&
	  firstState.action == BOBOL_ACTION_VIEW_PRIMARY_RELEASE,
	  "action-layer teardown is owner-conditional and restores view input");

    BObolInputBinding mged_bindings[] = {
	{BOBOL_INPUT_POINTER_PRESS, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_CONTROL, BOBOL_INPUT_MOD_SHIFT |
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_ROTATE_BEGIN},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_CONTROL, BOBOL_INPUT_MOD_SHIFT |
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_ROTATE_CANCEL},
	{BOBOL_INPUT_POINTER_PRESS, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_SHIFT, BOBOL_INPUT_MOD_CONTROL |
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_TRANSLATE_BEGIN},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_SHIFT, BOBOL_INPUT_MOD_CONTROL |
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_TRANSLATE_CANCEL},
	{BOBOL_INPUT_POINTER_PRESS, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_SHIFT | BOBOL_INPUT_MOD_CONTROL,
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_SCALE_BEGIN},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_SHIFT | BOBOL_INPUT_MOD_CONTROL,
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_SCALE_CANCEL},
	{BOBOL_INPUT_POINTER_PRESS, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY,
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_SHIFT,
	 BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_CONSTRAINED_TRANSLATE_BEGIN},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY,
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_SHIFT,
	 BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_CONSTRAINED_TRANSLATE_CANCEL},
	{BOBOL_INPUT_POINTER_PRESS, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY,
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_CONTROL,
	 BOBOL_INPUT_MOD_SHIFT | BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_CONSTRAINED_ROTATE_BEGIN},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY,
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_CONTROL,
	 BOBOL_INPUT_MOD_SHIFT | BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_CONSTRAINED_ROTATE_CANCEL},
	{BOBOL_INPUT_POINTER_PRESS, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY,
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_SHIFT |
	 BOBOL_INPUT_MOD_CONTROL, BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_CONSTRAINED_SCALE_BEGIN},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY,
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_SHIFT |
	 BOBOL_INPUT_MOD_CONTROL, BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_CONSTRAINED_SCALE_CANCEL},
	{BOBOL_INPUT_KEY_PRESS, BOBOL_INPUT_KEY_F3, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_SHIFT | BOBOL_INPUT_MOD_CONTROL |
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_PERSPECTIVE_MODE_TOGGLE},
	{BOBOL_INPUT_KEY_PRESS, BOBOL_INPUT_KEY_F6, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_SHIFT | BOBOL_INPUT_MOD_CONTROL |
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_PERSPECTIVE_CYCLE},
	{BOBOL_INPUT_KEY_PRESS, BOBOL_INPUT_KEY_F7, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_SHIFT | BOBOL_INPUT_MOD_CONTROL |
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_FACEPLATE_TOGGLE},
	{BOBOL_INPUT_KEY_PRESS, BOBOL_INPUT_KEY_F8, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_SHIFT | BOBOL_INPUT_MOD_CONTROL |
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META, 10,
	 TEST_LAYER_FACEPLATE_GUI_TOGGLE}
    };
    const BObolInputActionLayer mged_layer = {"test-mged-active-mode",
	mged_bindings, sizeof(mged_bindings) / sizeof(mged_bindings[0]),
	test_input_action};
    semantic_state = InputActionState();
    CHECK(first.setActionLayer(&mged_layer, &semantic_state, &semantic_state),
	  "application lifecycles install as an opaque action layer");
    event.type = BOBOL_INPUT_POINTER_PRESS;
    event.button = 0;
    event.modifiers = BOBOL_INPUT_MOD_CONTROL;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_HANDLED &&
	  semantic_state.action == TEST_LAYER_ROTATE_BEGIN,
	  "MGED rotate begin preserves a typed handled result");
    semantic_state.result = BOBOL_INPUT_RESULT_CANCELLED;
    event.type = BOBOL_INPUT_POINTER_RELEASE;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_CANCELLED &&
	  semantic_state.action == TEST_LAYER_ROTATE_CANCEL,
	  "MGED rotate release preserves a typed cancellation result");
    semantic_state.result = BOBOL_INPUT_RESULT_HANDLED;
    event.type = BOBOL_INPUT_POINTER_PRESS;
    event.modifiers = BOBOL_INPUT_MOD_SHIFT;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_HANDLED &&
	  semantic_state.action == TEST_LAYER_TRANSLATE_BEGIN,
	  "MGED translation begin preserves a typed handled result");
    semantic_state.result = BOBOL_INPUT_RESULT_CANCELLED;
    event.type = BOBOL_INPUT_POINTER_RELEASE;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_CANCELLED &&
	  semantic_state.action == TEST_LAYER_TRANSLATE_CANCEL,
	  "MGED translation release preserves a typed cancellation result");
    semantic_state.result = BOBOL_INPUT_RESULT_HANDLED;
    event.type = BOBOL_INPUT_POINTER_PRESS;
    event.modifiers = BOBOL_INPUT_MOD_SHIFT | BOBOL_INPUT_MOD_CONTROL;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_HANDLED &&
	  semantic_state.action == TEST_LAYER_SCALE_BEGIN,
	  "MGED scale begin preserves a typed handled result");
    semantic_state.result = BOBOL_INPUT_RESULT_CANCELLED;
    event.type = BOBOL_INPUT_POINTER_RELEASE;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_CANCELLED &&
	  semantic_state.action == TEST_LAYER_SCALE_CANCEL,
	  "MGED scale release preserves a typed cancellation result");
    semantic_state.result = BOBOL_INPUT_RESULT_HANDLED;
    event.type = BOBOL_INPUT_POINTER_PRESS;
    event.button = 1;
    event.modifiers = BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_SHIFT;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_HANDLED &&
	  semantic_state.action ==
	  TEST_LAYER_CONSTRAINED_TRANSLATE_BEGIN,
	  "MGED constrained translation begin is a typed action");
    semantic_state.result = BOBOL_INPUT_RESULT_CANCELLED;
    event.type = BOBOL_INPUT_POINTER_RELEASE;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_CANCELLED &&
	  semantic_state.action ==
	  TEST_LAYER_CONSTRAINED_TRANSLATE_CANCEL,
	  "MGED constrained translation cancel is a typed result");
    semantic_state.result = BOBOL_INPUT_RESULT_HANDLED;
    event.type = BOBOL_INPUT_POINTER_PRESS;
    event.modifiers = BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_CONTROL;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_HANDLED &&
	  semantic_state.action ==
	  TEST_LAYER_CONSTRAINED_ROTATE_BEGIN,
	  "MGED constrained rotation begin is a typed action");
    semantic_state.result = BOBOL_INPUT_RESULT_CANCELLED;
    event.type = BOBOL_INPUT_POINTER_RELEASE;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_CANCELLED &&
	  semantic_state.action ==
	  TEST_LAYER_CONSTRAINED_ROTATE_CANCEL,
	  "MGED constrained rotation cancel is a typed result");
    semantic_state.result = BOBOL_INPUT_RESULT_HANDLED;
    event.type = BOBOL_INPUT_POINTER_PRESS;
    event.modifiers = BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_SHIFT |
	BOBOL_INPUT_MOD_CONTROL;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_HANDLED &&
	  semantic_state.action ==
	  TEST_LAYER_CONSTRAINED_SCALE_BEGIN,
	  "MGED constrained scale begin is a typed action");
    semantic_state.result = BOBOL_INPUT_RESULT_CANCELLED;
    event.type = BOBOL_INPUT_POINTER_RELEASE;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_CANCELLED &&
	  semantic_state.action ==
	  TEST_LAYER_CONSTRAINED_SCALE_CANCEL,
	  "MGED constrained scale cancel is a typed result");
    semantic_state.result = BOBOL_INPUT_RESULT_HANDLED;
    event.type = BOBOL_INPUT_KEY_PRESS;
    event.key = BOBOL_INPUT_KEY_F3;
    event.button = BOBOL_INPUT_ANY;
    event.modifiers = BOBOL_INPUT_MOD_NONE;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_HANDLED &&
	  semantic_state.action == TEST_LAYER_PERSPECTIVE_MODE_TOGGLE,
	  "application projection-mode policy uses an opaque local action");
    event.key = BOBOL_INPUT_KEY_F6;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_HANDLED &&
	  semantic_state.action == TEST_LAYER_PERSPECTIVE_CYCLE,
	  "application projection-cycle policy uses an opaque local action");
    event.key = BOBOL_INPUT_KEY_F7;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_HANDLED &&
	  semantic_state.action == TEST_LAYER_FACEPLATE_TOGGLE,
	  "toolkit-neutral function keys dispatch application-local faceplate actions");
    event.key = BOBOL_INPUT_KEY_F8;
    CHECK(first.dispatch(&event) == BOBOL_INPUT_RESULT_HANDLED &&
	  semantic_state.action == TEST_LAYER_FACEPLATE_GUI_TOGGLE,
	  "adjacent application faceplate policy remains outside the standard action vocabulary");
    CHECK(first.clearActionLayerIf(&semantic_state),
	  "application lifecycles tear down by owner");

    BObolInputContext keyboard;
    keyboard.setProfile(bobol_input_keyboard_view_profile());
    keyboard.setActionHandler(test_input_action, &secondState);
    event.type = BOBOL_INPUT_KEY_PRESS;
    event.key = 'F';
    event.modifiers = BOBOL_INPUT_MOD_NONE;
    CHECK(keyboard.dispatch(&event) == 1 &&
	  secondState.action == BOBOL_ACTION_VIEW_FRONT,
	  "keyboard profile dispatches shared view shortcuts");
    event.type = BOBOL_INPUT_WHEEL;
    CHECK(keyboard.dispatch(&event) == 0,
	  "keyboard profile leaves application pointer gestures unhandled");

    struct bv_context *orientation_context = bv_context_create();
    vect_t aet = VINIT_ZERO;
    CHECK(orientation_context &&
	  bobol_input_view_orientation_apply(orientation_context,
	      BOBOL_ACTION_VIEW_TOP) &&
	  bv_aet_get(aet, bv_context_view(orientation_context)) &&
	  float_equal((float)aet[0], 270.0f) &&
	  float_equal((float)aet[1], 90.0f) &&
	  float_equal((float)aet[2], 0.0f),
	  "shared orientation helper applies the canonical top view");
    CHECK(!bobol_input_view_orientation_apply(orientation_context,
	      BOBOL_ACTION_VIEW_ZOOM),
	  "shared orientation helper leaves non-orientation actions unhandled");

    struct bv_adc_state adc = {};
    int visible = -1;
    CHECK(bv_adc_state_get(&adc, bv_context_view(orientation_context)) &&
	  bobol_display_endpoint_input_faceplate_toggle_apply(NULL,
	      orientation_context, BOBOL_ACTION_TOGGLE_ADC, &visible) &&
	  visible == !adc.draw &&
	  bv_adc_state_get(&adc, bv_context_view(orientation_context)) &&
	  adc.draw == visible,
	  "shared faceplate action toggles standalone ADC state");
    struct bv_axes_state model_axes = {};
    visible = -1;
    CHECK(bv_model_axes_state_get(&model_axes,
	  bv_context_view(orientation_context)) &&
	  bobol_display_endpoint_input_faceplate_toggle_apply(NULL,
	      orientation_context, BOBOL_ACTION_TOGGLE_MODEL_AXES, &visible) &&
	  visible == !model_axes.draw &&
	  bv_model_axes_state_get(&model_axes,
	      bv_context_view(orientation_context)) && model_axes.draw == visible,
	  "shared faceplate action toggles standalone model axes");
    struct bv_axes_state view_axes = {};
    visible = -1;
    CHECK(bv_view_axes_state_get(&view_axes,
	  bv_context_view(orientation_context)) &&
	  bobol_display_endpoint_input_faceplate_toggle_apply(NULL,
	      orientation_context, BOBOL_ACTION_TOGGLE_VIEW_AXES, &visible) &&
	  visible == !view_axes.draw &&
	  bv_view_axes_state_get(&view_axes,
	      bv_context_view(orientation_context)) && view_axes.draw == visible,
	  "shared faceplate action toggles standalone view axes");
    CHECK(!bobol_display_endpoint_input_faceplate_toggle_apply(NULL,
	orientation_context, BOBOL_ACTION_VIEW_TOP, NULL),
	"shared faceplate action rejects non-faceplate input");
    bv_context_destroy(orientation_context);
    return 0;
}

static int
test_progressive_status_contract(void)
{
    BObolViewController controller;
    PendingProgressState state;
    const uint64_t token = controller.registerProgressiveProvider(
	test_pending_progress_provider, &state);
    CHECK(token != 0, "progressive provider registers");

    BObolProgressiveStatus status;
    CHECK(controller.advanceProgressiveWork(NULL, &status) > 0 &&
	  state.calls == 1 && status.providerAdvanced == 1 &&
	  status.hasMore && !status.changed,
	  "pending progressive work does not claim a published scene change");
    controller.unregisterProgressiveProvider(token);
    controller.clearProgressiveWorkPending();
    return 0;
}

static int
test_window_host_contract(void)
{
    BObolWindowHost host;
    BObolWindowDesc desc;

    desc.mode = BOBOL_WINDOW_HEADLESS;
    desc.backend = BOBOL_WINDOW_BACKEND_OFFSCREEN;
    desc.width = 6;
    desc.height = 4;
    desc.title = "window-host-test";
    desc.display = "";
    desc.nativeIdHint = "";
    desc.visible = FALSE;

    CHECK(host.open(&desc) == 0, "window host opens neutral controller");
    CHECK(host.isOpen(), "window host records open state");
    CHECK(host.getController() != NULL, "window host owns a view controller");
    CHECK(host.getController()->getSceneRoot() != NULL, "window host creates scene root");
    CHECK(host.getController()->getSceneRoot()->isOfType(SoGroup::getClassTypeId()),
	  "window host scene root is a group");
    CHECK(host.getController()->getViewportRegion().getWindowSize()[0] == 6 &&
	  host.getController()->getViewportRegion().getWindowSize()[1] == 4,
	  "window host applies requested viewport size");

    BObolInputBinding binding;
    binding.eventType = BOBOL_INPUT_KEY_PRESS;
    binding.key = 'f';
    binding.button = BOBOL_INPUT_ANY;
    binding.requiredModifiers = 0;
    binding.action = BOBOL_ACTION_VIEW_FRONT;

    BObolInputProfile profile;
    profile.name = "test";
    profile.bindings = &binding;
    profile.bindingCount = 1;

    BObolInputEvent event;
    event.type = BOBOL_INPUT_KEY_PRESS;
    event.key = 'f';
    CHECK(host.handleInputEvent(&event, &profile) == 1,
	  "window host applies semantic input profile action");
    CHECK(host.getController()->isRenderRequested(),
	  "semantic input action requests render");
    host.getController()->clearRenderRequest();

    host.close();
    CHECK(!host.isOpen(), "window host closes");
    return 0;
}

class CountingWindowHost : public BObolWindowHost {
public:
    ~CountingWindowHost(void) override
    {
	destroyed++;
    }

    static int destroyed;
};

int CountingWindowHost::destroyed = 0;

struct FactoryTestState {
    FactoryTestState(void) :
	probe_result(1), open_result(1), creates(0), destroys(0), binds(0),
	detaches(0), opens(0), closes(0), frames(0), resizes(0), visible(0),
	vsync(-1), block_frames(false), frame_entered(false),
	frame_release(false), frame_returned(false)
    {
    }

    int probe_result;
    int open_result;
    int creates;
    int destroys;
    int binds;
    int detaches;
    int opens;
    int closes;
    int frames;
    int resizes;
    int captures = 0;
    int dimension_queries = 0;
    int framebuffer_captures = 0;
    std::string title;
    int visible;
    int vsync;
    BObolInputEventHandler input_dispatch = NULL;
    void *input_dispatch_data = NULL;
    std::mutex frame_mutex;
    std::condition_variable frame_cv;
    bool block_frames;
    bool frame_entered;
    bool frame_release;
    bool frame_returned;
};

struct FactoryTestInstance {
    FactoryTestState *state;
    void *controller;
    BObolInputEventHandler input_dispatch;
    void *input_dispatch_data;
    BObolWindowHost framebuffer_host;
};

static int
factory_test_probe(const struct bobol_host_desc *UNUSED(desc), void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    return state->probe_result;
}

static void *
factory_test_create(const struct bobol_host_desc *desc, void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    FactoryTestInstance *instance = new FactoryTestInstance;
    instance->state = state;
    instance->controller = NULL;
    instance->input_dispatch = desc ? desc->input_dispatch : NULL;
    instance->input_dispatch_data = desc ? desc->input_dispatch_data : NULL;
    state->input_dispatch = desc ? desc->input_dispatch : NULL;
    state->input_dispatch_data = desc ? desc->input_dispatch_data : NULL;
    state->creates++;
    return instance;
}

static void
factory_test_destroy(void *instance_ptr, void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    state->destroys++;
    delete static_cast<FactoryTestInstance *>(instance_ptr);
}

static int
factory_test_bind(void *instance_ptr, void *controller, void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    FactoryTestInstance *instance =
	static_cast<FactoryTestInstance *>(instance_ptr);
    BObolViewController *oldController =
	static_cast<BObolViewController *>(instance->controller);
    BObolViewController *newController =
	static_cast<BObolViewController *>(controller);
    SoDB::ContextManager *provider = bobol_headless_context_manager();
    if (oldController && oldController != newController &&
	oldController->getRenderContextManager() == provider)
	oldController->setRenderContextManager(NULL);
    if (oldController != newController) {
	instance->controller = controller;
	instance->framebuffer_host.attachController(newController, FALSE);
	if (controller)
	    state->binds++;
	else
	    state->detaches++;
    }
    if (newController)
	newController->setRenderContextManager(provider);
    return 1;
}

static int
factory_test_bind_without_provider(void *instance_ptr, void *controller,
	void *data)
{
    const int ret = factory_test_bind(instance_ptr, controller, data);
    BObolViewController *bound =
	static_cast<BObolViewController *>(controller);
    if (bound)
	bound->setRenderContextManager(NULL);
    return ret;
}

static void *
factory_test_framebuffer_window_host(void *instance_ptr,
	void *UNUSED(data))
{
    FactoryTestInstance *instance =
	static_cast<FactoryTestInstance *>(instance_ptr);
    return instance ? &instance->framebuffer_host : NULL;
}

static int
factory_test_open(void *UNUSED(instance),
	const struct bobol_host_desc *UNUSED(desc), void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    state->opens++;
    return state->open_result;
}

static void
factory_test_close(void *UNUSED(instance), void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    state->closes++;
}

static int
factory_test_frame(void *UNUSED(instance), const char *UNUSED(reason),
	void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    std::unique_lock<std::mutex> lock(state->frame_mutex);
    state->frames++;
    if (state->block_frames) {
	state->frame_entered = true;
	state->frame_cv.notify_all();
	state->frame_cv.wait(lock, [state]() {
	    return state->frame_release;
	});
	state->frame_returned = true;
	state->frame_cv.notify_all();
    }
    return 1;
}

static int
factory_test_resize(void *UNUSED(instance), unsigned int UNUSED(width),
	unsigned int UNUSED(height), double UNUSED(device_pixel_ratio), void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    state->resizes++;
    return 1;
}

static int
factory_test_capture(void *UNUSED(instance), unsigned char **pixels,
	size_t *size, unsigned int *width, unsigned int *height,
	unsigned int *components, void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    *width = 2;
    *height = 2;
    *components = 3;
    *size = 12;
    *pixels = static_cast<unsigned char *>(bu_calloc(*size, 1,
	"factory test capture"));
    state->captures++;
    return 1;
}

static int
factory_test_dimensions(void *UNUSED(instance), unsigned int *width,
	unsigned int *height, double *device_pixel_ratio, void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    *width = 13;
    *height = 9;
    *device_pixel_ratio = 2.0;
    state->dimension_queries++;
    return 1;
}

static int
factory_test_set_title(void *UNUSED(instance), const char *title, void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    if (!title)
	return 0;
    state->title = title;
    return 1;
}

static int
factory_test_set_visible(void *UNUSED(instance), int visible, void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    state->visible = visible ? 1 : 0;
    return 1;
}

static int
factory_test_set_vsync(void *UNUSED(instance), int enabled, void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    state->vsync = enabled ? 1 : 0;
    return 1;
}

static int
factory_test_framebuffer_capture(void *data, unsigned char **pixels,
	size_t *size, unsigned int *width, unsigned int *height,
	unsigned int *components)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    if (!state || !pixels || !size || !width || !height || !components)
	return 0;
    *width = 1;
    *height = 1;
    *components = 3;
    *size = 3;
    *pixels = static_cast<unsigned char *>(bu_malloc(*size,
	"factory framebuffer capture"));
    (*pixels)[0] = 17;
    (*pixels)[1] = 34;
    (*pixels)[2] = 51;
    state->framebuffer_captures++;
    return 1;
}

static struct bobol_host_factory
factory_test_desc(const char *name, int priority, uint64_t capabilities,
	FactoryTestState *state)
{
    struct bobol_host_factory factory;
    memset(&factory, 0, sizeof(factory));
    factory.abi_version = BOBOL_HOST_FACTORY_ABI_VERSION;
    factory.struct_size = sizeof(factory);
    factory.name = name;
    factory.priority = priority;
    factory.capabilities = capabilities;
    factory.user_data = state;
    factory.probe = factory_test_probe;
    factory.create = factory_test_create;
    factory.destroy = factory_test_destroy;
    factory.bind_controller = factory_test_bind;
    factory.open = factory_test_open;
    factory.close = factory_test_close;
    factory.request_frame = factory_test_frame;
    factory.resize = factory_test_resize;
    factory.capture = factory_test_capture;
    factory.dimensions = factory_test_dimensions;
    factory.set_title = factory_test_set_title;
    factory.set_visible = factory_test_set_visible;
    factory.set_vsync = factory_test_set_vsync;
    return factory;
}

static int
test_host_factory_contract(void)
{
    const size_t initial_count = bobol_host_factory_registry_count();
    FactoryTestState low_state;
    FactoryTestState high_state;
    FactoryTestState failed_state;
    FactoryTestState missing_provider_state;
    failed_state.open_result = 0;

    struct bobol_host_factory invalid =
	factory_test_desc("endpoint-test-invalid", 0, 0, &low_state);
    invalid.abi_version++;
    CHECK(bobol_host_factory_register(&invalid) == NULL,
	  "host factory rejects an unsupported ABI");

    struct bobol_host_factory low = factory_test_desc(
	"endpoint-test-low", 10,
	BOBOL_HOST_CAP_SYSTEM_GL | BOBOL_HOST_CAP_PIXEL_PRESENT,
	&low_state);
    struct bobol_host_factory high = factory_test_desc(
	"endpoint-test-high", 20,
	BOBOL_HOST_CAP_TOPLEVEL | BOBOL_HOST_CAP_PIXEL_PRESENT |
	BOBOL_HOST_CAP_READBACK | BOBOL_HOST_CAP_PRESENT_VSYNC,
	&high_state);
    high.framebuffer_window_host = factory_test_framebuffer_window_host;
    struct bobol_host_factory failed = factory_test_desc(
	"endpoint-test-failed", 30, BOBOL_HOST_CAP_PIXEL_PRESENT,
	&failed_state);
    struct bobol_host_factory missing_provider = factory_test_desc(
	"endpoint-test-missing-provider", 30, BOBOL_HOST_CAP_PIXEL_PRESENT,
	&missing_provider_state);
    missing_provider.bind_controller = factory_test_bind_without_provider;

    bobol_host_factory_token_t *low_token =
	bobol_host_factory_register(&low);
    bobol_host_factory_token_t *high_token =
	bobol_host_factory_register(&high);
    bobol_host_factory_token_t *failed_token =
	bobol_host_factory_register(&failed);
    bobol_host_factory_token_t *missing_provider_token =
	bobol_host_factory_register(&missing_provider);
    CHECK(low_token && high_token && failed_token && missing_provider_token,
	  "host factories register copied descriptors");
    CHECK(bobol_host_factory_register(&high) == NULL,
	  "host factory rejects duplicate stable names");
    CHECK(bobol_host_factory_registry_count() == initial_count + 4,
	  "host factory registry reports registrations");
    int found_high_capabilities = 0;
    for (size_t i = 0; i < bobol_host_factory_registry_count(); i++) {
	char name[64] = {0};
	if (bobol_host_factory_registry_name(i, name, sizeof(name)) &&
	    bu_strcmp(name, "endpoint-test-high") == 0) {
	    found_high_capabilities =
		bobol_host_factory_registry_capabilities(i) ==
		(BOBOL_HOST_CAP_PIXEL_PRESENT |
		 BOBOL_HOST_CAP_READBACK | BOBOL_HOST_CAP_TOPLEVEL |
		 BOBOL_HOST_CAP_PRESENT_VSYNC);
	}
    }
    CHECK(found_high_capabilities,
	  "host factory registry exposes immutable capabilities");

    struct bobol_host_desc desc;
    memset(&desc, 0, sizeof(desc));
    desc.struct_size = sizeof(desc);
    desc.mode = BOBOL_HOST_MODE_TOPLEVEL;
    desc.width = 8;
    desc.height = 6;
    desc.device_pixel_ratio = 1.0;
    desc.vsync = BOBOL_HOST_VSYNC_OFF;
    desc.required_capabilities = BOBOL_HOST_CAP_READBACK;

    bobol_display_endpoint_t *endpoint =
	bobol_display_endpoint_create(NULL, 0);
    desc.required_capabilities = BOBOL_HOST_CAP_PIXEL_PRESENT;
    desc.vsync = BOBOL_HOST_VSYNC_AUTO;
    CHECK(endpoint && !bobol_display_endpoint_host_open(endpoint,
	  "endpoint-test-missing-provider", &desc) &&
	  missing_provider_state.opens == 1 &&
	  missing_provider_state.closes == 1 &&
	  missing_provider_state.detaches == 1 &&
	  missing_provider_state.destroys == 1,
	  "endpoint rejects a rendering factory that binds no provider");
    bobol_display_endpoint_destroy(endpoint);

    desc.required_capabilities = BOBOL_HOST_CAP_READBACK;
    desc.vsync = BOBOL_HOST_VSYNC_OFF;

    endpoint = bobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint != NULL, "factory test creates display endpoint");
    CHECK(bobol_display_endpoint_host_open(endpoint, NULL, &desc),
	  "endpoint selects a compatible registered factory");
    CHECK(bu_strcmp(bobol_display_endpoint_host_factory_name(endpoint),
	  "endpoint-test-high") == 0,
	  "endpoint selection honors capabilities before priority");
    CHECK(bobol_display_endpoint_host_capabilities(endpoint) ==
	  (BOBOL_HOST_CAP_TOPLEVEL | BOBOL_HOST_CAP_PIXEL_PRESENT |
	   BOBOL_HOST_CAP_READBACK | BOBOL_HOST_CAP_PRESENT_VSYNC),
	  "endpoint exposes its active host capabilities");
    CHECK(bobol_display_endpoint_render_engine_capabilities(endpoint) ==
	  (BOBOL_RENDER_ENGINE_CAP_AUTO | BOBOL_RENDER_ENGINE_CAP_SW |
	   BOBOL_RENDER_ENGINE_CAP_NONE |
	   BOBOL_RENDER_ENGINE_CAP_DIAGNOSTIC) &&
	  bobol_display_endpoint_render_engine_supported(endpoint,
	      BOBOL_RENDER_ENGINE_SW) &&
	  !bobol_display_endpoint_render_engine_supported(endpoint,
	      BOBOL_RENDER_ENGINE_HW),
	  "endpoint reports renderer support from its active host capabilities");
    struct bv_display_property_value resolved_property =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    CHECK(bobol_display_endpoint_render_engine_get(endpoint) ==
	  BOBOL_RENDER_ENGINE_AUTO &&
	  bobol_display_endpoint_render_engine_resolved_get(endpoint) ==
	  BOBOL_RENDER_ENGINE_SW &&
	  bobol_display_endpoint_property_get(endpoint,
	      "endpoint.renderer.resolved", &resolved_property) ==
	      BV_DISPLAY_PROPERTY_OK &&
	  resolved_property.string_value &&
	  bu_strcmp(resolved_property.string_value, "sw") == 0,
	  "automatic renderer resolves a pixel-presentation factory to software");
    CHECK(!bobol_display_endpoint_render_engine_set(endpoint,
	  BOBOL_RENDER_ENGINE_RT),
	  "retained rt rejects a host without progressive presentation");
    FactoryTestInstance *instance = static_cast<FactoryTestInstance *>(
	bobol_display_endpoint_host(endpoint));
    CHECK(instance && instance->controller ==
	  bobol_display_endpoint_controller(endpoint),
	  "factory host instance binds the endpoint controller");
    CHECK(bobol_display_endpoint_framebuffer_window_host(endpoint) ==
	  &instance->framebuffer_host &&
	  instance->framebuffer_host.getController() ==
	  bobol_display_endpoint_controller(endpoint),
	  "factory exposes its retained framebuffer host through the endpoint");
    InputActionState input_state;
    BObolInputEvent input_event;
    input_event.type = BOBOL_INPUT_KEY_PRESS;
    input_event.key = 'F';
    CHECK(bobol_display_endpoint_input_action_handler_set(endpoint,
	  test_input_action, &input_state) && high_state.input_dispatch &&
	  high_state.input_dispatch_data == endpoint &&
	  high_state.input_dispatch(high_state.input_dispatch_data, &input_event) == 1 &&
	  input_state.action == BOBOL_ACTION_VIEW_FRONT,
	  "factory receives the endpoint-local semantic input dispatcher");
    CHECK(!bobol_host_factory_unregister(high_token),
	  "live endpoint prevents host factory unregister");

    struct bv_display_property_value host_dimension =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    CHECK(bobol_display_endpoint_property_get(endpoint, "endpoint.width",
	  &host_dimension) == BV_DISPLAY_PROPERTY_OK &&
	  host_dimension.uint_value == 13 && high_state.dimension_queries == 1,
	  "endpoint dimensions refresh from the active toolkit host");

    struct bv_display_property_value host_property =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    host_property.type = BV_DISPLAY_PROPERTY_STRING;
    host_property.string_value = "Updated endpoint title";
    CHECK(bobol_display_endpoint_property_set(endpoint, "endpoint.title",
	  &host_property) == BV_DISPLAY_PROPERTY_OK &&
	  high_state.title == "Updated endpoint title",
	  "typed title property dispatches to a toplevel host");
    property_value_init(&host_property);
    CHECK(bobol_display_endpoint_property_get(endpoint, "endpoint.title",
	  &host_property) == BV_DISPLAY_PROPERTY_OK &&
	  host_property.string_value &&
	  bu_strcmp(host_property.string_value, "Updated endpoint title") == 0,
	  "typed title property retains endpoint state");
    property_value_init(&host_property);
    host_property.type = BV_DISPLAY_PROPERTY_BOOL;
    host_property.bool_value = 1;
    CHECK(bobol_display_endpoint_property_set(endpoint, "endpoint.visible",
	  &host_property) == BV_DISPLAY_PROPERTY_OK && high_state.visible,
	  "typed visibility property dispatches to a toplevel host");
    property_value_init(&host_property);
    CHECK(bobol_display_endpoint_property_get(endpoint, "endpoint.vsync",
	  &host_property) == BV_DISPLAY_PROPERTY_OK &&
	  host_property.bool_value == 0,
	  "explicit host creation policy initializes typed vsync state");
    host_property.type = BV_DISPLAY_PROPERTY_BOOL;
    host_property.bool_value = 1;
    CHECK(bobol_display_endpoint_property_set(endpoint, "endpoint.vsync",
	  &host_property) == BV_DISPLAY_PROPERTY_OK &&
	  high_state.vsync == 1,
	  "typed vsync property dispatches to a capable presentation host");

	const int frames_before_request = high_state.frames;
    CHECK(bobol_display_endpoint_request_frame(endpoint, "test-frame") &&
	  high_state.frames == frames_before_request + 1,
	  "factory dispatches frame requests");
    CHECK(bobol_display_endpoint_resize(endpoint, 10, 7, 1.5) &&
	  high_state.resizes == 1 &&
	  instance->controller == bobol_display_endpoint_controller(endpoint),
	  "factory dispatches resize requests");
    unsigned char *capture = NULL;
    size_t capture_size = 0;
    unsigned int capture_width = 0;
    unsigned int capture_height = 0;
    unsigned int capture_components = 0;
    CHECK(bobol_display_endpoint_capture(endpoint, &capture, &capture_size,
	  &capture_width, &capture_height, &capture_components) &&
	  capture && capture_size == 12 && capture_width == 2 &&
	  capture_height == 2 && capture_components == 3 &&
	  high_state.captures == 1,
	  "endpoint dispatches capture through a readback factory");
    bu_free(capture, "factory test capture");

    CHECK(!bobol_display_endpoint_capture_plane(endpoint,
	  BOBOL_CAPTURE_FRAMEBUFFER, &capture, &capture_size,
	  &capture_width, &capture_height, &capture_components),
	  "framebuffer-plane capture fails without an explicit provider");
    CHECK(bobol_display_endpoint_framebuffer_capture_provider_set(endpoint,
	  factory_test_framebuffer_capture, &high_state),
	  "endpoint accepts a per-instance framebuffer capture provider");
    CHECK(bobol_display_endpoint_capture_plane(endpoint,
	  BOBOL_CAPTURE_FRAMEBUFFER, &capture, &capture_size,
	  &capture_width, &capture_height, &capture_components) && capture &&
	  capture_size == 3 && capture_width == 1 && capture_height == 1 &&
	  capture_components == 3 && capture[0] == 17 && capture[1] == 34 &&
	  capture[2] == 51 && high_state.framebuffer_captures == 1,
	  "endpoint dispatches framebuffer-plane capture to its provider");
    bu_free(capture, "factory framebuffer capture");

    BObolViewController *host_controller =
	static_cast<BObolViewController *>(
	    bobol_display_endpoint_controller(endpoint));
    SoNode *graphical_root = host_controller->getRenderRoot();
    const int frames_before_none = high_state.frames;
    const int captures_before_none = high_state.captures;
    CHECK(graphical_root && bobol_display_endpoint_render_engine_set(
	  endpoint, BOBOL_RENDER_ENGINE_NONE) &&
	  host_controller->getSceneRoot() && !host_controller->getRenderRoot() &&
	  bobol_display_endpoint_request_frame(endpoint, "none-test") &&
	  high_state.frames == frames_before_none &&
	  !bobol_display_endpoint_capture(endpoint, &capture, &capture_size,
	      &capture_width, &capture_height, &capture_components) &&
	  high_state.captures == captures_before_none,
	  "none retains the scene while suppressing hosted frame and capture paths");
    CHECK(bobol_display_endpoint_capture_plane(endpoint,
	  BOBOL_CAPTURE_FRAMEBUFFER, &capture, &capture_size,
	  &capture_width, &capture_height, &capture_components) && capture &&
	  capture[0] == 17 && capture[1] == 34 && capture[2] == 51,
	  "none preserves the independent framebuffer capture plane");
    bu_free(capture, "none framebuffer capture");

    CHECK(bobol_display_endpoint_render_engine_set(endpoint,
	  BOBOL_RENDER_ENGINE_DIAGNOSTIC) &&
	  !host_controller->getRenderRoot(),
	  "diagnostic disables graphical presentation on a hosted endpoint");
    struct bv_display_property_value diagnostic_value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    CHECK(bobol_display_endpoint_property_get(endpoint,
	  "renderer.diagnostic.revision", &diagnostic_value) ==
	  BV_DISPLAY_PROPERTY_OK && diagnostic_value.uint_value > 0,
	  "diagnostic selection performs a typed non-graphical traversal");
    const uint64_t diagnostic_revision = diagnostic_value.uint_value;
    property_value_init(&diagnostic_value);
    CHECK(bobol_display_endpoint_request_frame(endpoint, "diagnostic-test") &&
	  high_state.frames == frames_before_none &&
	  bobol_display_endpoint_property_get(endpoint,
	      "renderer.diagnostic.revision", &diagnostic_value) ==
	      BV_DISPLAY_PROPERTY_OK &&
	  diagnostic_value.uint_value > diagnostic_revision &&
	  !bobol_display_endpoint_capture(endpoint, &capture, &capture_size,
	      &capture_width, &capture_height, &capture_components) &&
	  high_state.captures == captures_before_none,
	  "diagnostic refreshes reports without posing as a graphical backend");
    property_value_init(&diagnostic_value);
    CHECK(bobol_display_endpoint_property_get(endpoint,
	  "renderer.diagnostic.summary", &diagnostic_value) ==
	  BV_DISPLAY_PROPERTY_OK && diagnostic_value.string_value &&
	  strstr(diagnostic_value.string_value, "scene=available") &&
	  strstr(diagnostic_value.string_value, "visited_sources="),
	  "diagnostic exposes its traversal result through a typed summary");
    CHECK(bobol_display_endpoint_render_engine_set(endpoint,
	  BOBOL_RENDER_ENGINE_SW) &&
	  bobol_display_endpoint_render_engine_resolved_get(endpoint) ==
	  BOBOL_RENDER_ENGINE_SW &&
	  host_controller->getRenderRoot() == graphical_root &&
	  bobol_display_endpoint_property_get(endpoint,
	      "renderer.diagnostic.summary", &diagnostic_value) ==
	      BV_DISPLAY_PROPERTY_UNSUPPORTED,
	  "leaving diagnostic restores the retained graphical root and hides mode-only state");

    CHECK(!bobol_display_endpoint_framebuffer_capture_provider_set(endpoint,
	  NULL, &low_state),
	  "unrelated owner cannot clear an endpoint framebuffer provider");
    CHECK(bobol_display_endpoint_framebuffer_capture_provider_set(endpoint,
	  NULL, &high_state),
	  "framebuffer provider owner can clear its endpoint binding");

    bobol_display_endpoint_destroy(endpoint);
    CHECK(high_state.binds == 1 && high_state.opens == 1 &&
	  high_state.closes == 1 && high_state.detaches == 1 &&
	  high_state.destroys == 1,
	  "endpoint closes, detaches, and destroys its factory host");

    desc.required_capabilities = BOBOL_HOST_CAP_PIXEL_PRESENT;
    desc.vsync = BOBOL_HOST_VSYNC_AUTO;
    endpoint = bobol_display_endpoint_create(NULL, 0);
    CHECK(!bobol_display_endpoint_host_open(endpoint,
	  "endpoint-test-failed", &desc),
	  "endpoint reports explicit factory open failure");
    CHECK(failed_state.creates == 1 && failed_state.opens == 1 &&
	  failed_state.detaches == 1 && failed_state.destroys == 1 &&
	  failed_state.closes == 0,
	  "failed factory open detaches and destroys without close");
    bobol_display_endpoint_destroy(endpoint);

    desc.mode = BOBOL_HOST_MODE_HEADLESS;
    endpoint = bobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint && bobol_display_endpoint_host_open(endpoint,
	  "endpoint-test-low", &desc) &&
	  bobol_display_endpoint_render_engine_get(endpoint) ==
	  BOBOL_RENDER_ENGINE_AUTO &&
	  bobol_display_endpoint_render_engine_resolved_get(endpoint) ==
	  BOBOL_RENDER_ENGINE_HW,
	  "automatic renderer deterministically prefers System GL on a hybrid host");
    property_value_init(&resolved_property);
    CHECK(bobol_display_endpoint_property_get(endpoint,
	  "endpoint.renderer.resolved", &resolved_property) ==
	  BV_DISPLAY_PROPERTY_OK && resolved_property.string_value &&
	  bu_strcmp(resolved_property.string_value, "hw") == 0,
	  "typed resolved-renderer property reports automatic hardware selection");
    bobol_display_endpoint_destroy(endpoint);

    CHECK(bobol_host_factory_unregister(missing_provider_token) &&
	  bobol_host_factory_unregister(failed_token) &&
	  bobol_host_factory_unregister(high_token) &&
	  bobol_host_factory_unregister(low_token),
	  "unused host factories unregister");
    CHECK(bobol_host_factory_registry_count() == initial_count,
	  "host factory test restores registry state");
    return 0;
}

static int
test_host_factory_simultaneous_endpoints(void)
{
    FactoryTestState state;
    struct bobol_host_factory factory = factory_test_desc(
	"endpoint-test-simultaneous", 10,
	BOBOL_HOST_CAP_PIXEL_PRESENT | BOBOL_HOST_CAP_INPUT, &state);
    bobol_host_factory_token_t *token =
	bobol_host_factory_register(&factory);
    CHECK(token != NULL, "simultaneous endpoint factory registers");

    struct bobol_host_desc desc = {};
    desc.struct_size = sizeof(desc);
    desc.mode = BOBOL_HOST_MODE_HEADLESS;
    desc.width = 16;
    desc.height = 12;
    desc.device_pixel_ratio = 1.0;
    desc.required_capabilities = BOBOL_HOST_CAP_PIXEL_PRESENT;

    bobol_display_endpoint_t *first =
	bobol_display_endpoint_create(NULL, 0);
    bobol_display_endpoint_t *second =
	bobol_display_endpoint_create(NULL, 0);
    CHECK(first && second && bobol_display_endpoint_host_open(first,
	  "endpoint-test-simultaneous", &desc) &&
	  bobol_display_endpoint_host_open(second,
	  "endpoint-test-simultaneous", &desc),
	  "simultaneous endpoints open isolated factory instances");

    FactoryTestInstance *first_instance =
	static_cast<FactoryTestInstance *>(bobol_display_endpoint_host(first));
    FactoryTestInstance *second_instance =
	static_cast<FactoryTestInstance *>(bobol_display_endpoint_host(second));
    CHECK(first_instance && second_instance &&
	  first_instance != second_instance &&
	  first_instance->controller != second_instance->controller &&
	  first_instance->input_dispatch && second_instance->input_dispatch,
	  "simultaneous endpoints retain separate controllers and dispatchers");

    InputActionState first_input;
    InputActionState second_input;
    BObolInputEvent event;
    event.type = BOBOL_INPUT_KEY_PRESS;
    event.key = 'F';
    CHECK(bobol_display_endpoint_input_action_handler_set(first,
	  test_input_action, &first_input) &&
	  bobol_display_endpoint_input_action_handler_set(second,
	  test_input_action, &second_input) &&
	  first_instance->input_dispatch(first_instance->input_dispatch_data,
	  &event) == 1 && first_input.calls == 1 && second_input.calls == 0 &&
	  second_instance->input_dispatch(second_instance->input_dispatch_data,
	  &event) == 1 && first_input.calls == 1 && second_input.calls == 1,
	  "simultaneous endpoints dispatch input only to their own handler");

    BObolViewController *second_controller =
	static_cast<BObolViewController *>(
	bobol_display_endpoint_controller(second));
    const int frames_before = state.frames;
    CHECK(bobol_display_endpoint_request_frame(first, "first-frame") &&
	  bobol_display_endpoint_request_frame(second, "second-frame") &&
	  state.frames == frames_before + 2,
	  "simultaneous endpoints route their own frame requests");

    bobol_display_endpoint_destroy(first);
    const int frames_after_first_teardown = state.frames;
    second_controller->markProgressiveWorkPending();
    CHECK(state.frames == frames_after_first_teardown + 1 &&
	  state.closes == 1 && state.destroys == 1,
	  "surviving endpoint retains its progressive callback after peer teardown");
    second_controller->clearProgressiveWorkPending();

    bobol_display_endpoint_destroy(second);
    CHECK(state.closes == 2 && state.detaches == 2 && state.destroys == 2 &&
	  bobol_host_factory_unregister(token),
	  "simultaneous endpoint teardown releases every host instance");
    return 0;
}

static int
test_host_factory_reattach(void)
{
    FactoryTestState first_state;
    FactoryTestState second_state;
    second_state.open_result = 0;
    struct bobol_host_factory first_factory = factory_test_desc(
	"endpoint-test-reattach-first", 10,
	BOBOL_HOST_CAP_PIXEL_PRESENT | BOBOL_HOST_CAP_INPUT,
	&first_state);
    struct bobol_host_factory second_factory = factory_test_desc(
	"endpoint-test-reattach-second", 10,
	BOBOL_HOST_CAP_PIXEL_PRESENT | BOBOL_HOST_CAP_INPUT,
	&second_state);
    bobol_host_factory_token_t *first_token =
	bobol_host_factory_register(&first_factory);
    bobol_host_factory_token_t *second_token =
	bobol_host_factory_register(&second_factory);
    CHECK(first_token && second_token,
	  "reattach endpoint factories register");

    struct bobol_host_desc desc = {};
    desc.struct_size = sizeof(desc);
    desc.mode = BOBOL_HOST_MODE_HEADLESS;
    desc.width = 16;
    desc.height = 12;
    desc.device_pixel_ratio = 1.0;
    desc.required_capabilities = BOBOL_HOST_CAP_PIXEL_PRESENT;

    bobol_display_endpoint_t *endpoint =
	bobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint && bobol_display_endpoint_host_open(endpoint,
	  "endpoint-test-reattach-first", &desc),
	  "endpoint opens its initial factory host");
    FactoryTestInstance *first_instance =
	static_cast<FactoryTestInstance *>(
	bobol_display_endpoint_host(endpoint));
    void *controller = bobol_display_endpoint_controller(endpoint);
    CHECK(first_instance && controller && first_instance->controller == controller,
	  "initial factory host receives the endpoint controller");

    InputActionState input_state;
    BObolInputEvent input_event;
    input_event.type = BOBOL_INPUT_KEY_PRESS;
    input_event.key = 'F';
    CHECK(bobol_display_endpoint_input_action_handler_set(endpoint,
	  test_input_action, &input_state) && first_instance->input_dispatch &&
	  first_instance->input_dispatch(first_instance->input_dispatch_data,
	  &input_event) == 1 && input_state.calls == 1,
	  "initial factory host dispatches the endpoint-local input handler");

    const int first_frames = first_state.frames;
    CHECK(bobol_display_endpoint_request_frame(endpoint, "reattach-first") &&
	  first_state.frames == first_frames + 1,
	  "initial factory host receives endpoint frame requests");

    CHECK(!bobol_display_endpoint_host_open(endpoint,
	  "endpoint-test-reattach-second", &desc),
	  "failed replacement host leaves the existing endpoint host intact");
    CHECK(bobol_display_endpoint_host(endpoint) == first_instance &&
	  bu_strcmp(bobol_display_endpoint_host_factory_name(endpoint),
	  "endpoint-test-reattach-first") == 0 &&
	  static_cast<BObolViewController *>(controller)->
	      getRenderContextManager() == bobol_headless_context_manager() &&
	  first_state.closes == 0 && first_state.detaches == 0 &&
	  first_state.destroys == 0 && second_state.creates == 1 &&
	  second_state.binds == 1 && second_state.opens == 1 &&
	  second_state.closes == 0 && second_state.detaches == 1 &&
	  second_state.destroys == 1,
	  "failed replacement cleans only its provisional factory instance");
    CHECK(first_instance->input_dispatch(first_instance->input_dispatch_data,
	  &input_event) == 1 && input_state.calls == 2,
	  "failed replacement preserves the active host input route");

    second_state.open_result = 1;
    desc.width = 31;
    desc.height = 23;
    desc.device_pixel_ratio = 1.25;
    CHECK(bobol_display_endpoint_host_open(endpoint,
	  "endpoint-test-reattach-second", &desc),
	  "endpoint replaces its host after the new factory is ready");
    FactoryTestInstance *second_instance =
	static_cast<FactoryTestInstance *>(
	bobol_display_endpoint_host(endpoint));
    CHECK(second_instance && second_instance->controller == controller &&
	  bu_strcmp(bobol_display_endpoint_host_factory_name(endpoint),
	  "endpoint-test-reattach-second") == 0 &&
	  first_state.closes == 1 && first_state.detaches == 1 &&
	  first_state.destroys == 1 && second_state.creates == 2 &&
	  second_state.binds == 2 && second_state.opens == 2,
	  "replacement preserves endpoint identity and releases the old host");
    CHECK(second_instance->input_dispatch &&
	  second_instance->input_dispatch(second_instance->input_dispatch_data,
	  &input_event) == 1 && input_state.calls == 3,
	  "replacement host receives the existing endpoint input handler");

    const int first_frames_after_replace = first_state.frames;
    const int second_frames = second_state.frames;
    CHECK(bobol_display_endpoint_request_frame(endpoint, "reattach-second") &&
	  first_state.frames == first_frames_after_replace &&
	  second_state.frames == second_frames + 1,
	  "replacement routes future frame requests only to the active host");

    bobol_display_endpoint_host_detach(endpoint);
    CHECK(bobol_display_endpoint_host(endpoint) == NULL &&
	  second_state.closes == 1 && second_state.detaches == 2 &&
	  second_state.destroys == 2,
	  "endpoint detach releases the replacement host exactly once");
    bobol_display_endpoint_destroy(endpoint);
    CHECK(bobol_host_factory_unregister(second_token) &&
	  bobol_host_factory_unregister(first_token),
	  "reattach factories unregister after endpoint teardown");
    return 0;
}

struct BlockingFrameRequest {
    std::mutex mutex;
    std::condition_variable cv;
    bool entered = false;
    bool release = false;
    bool returned = false;
};

static void
blocking_frame_request(void *data, const char *UNUSED(reason))
{
    BlockingFrameRequest *state = static_cast<BlockingFrameRequest *>(data);
    std::unique_lock<std::mutex> lock(state->mutex);
    state->entered = true;
    state->cv.notify_all();
    state->cv.wait(lock, [state]() {
	return state->release;
    });
    state->returned = true;
    state->cv.notify_all();
}

static int
test_frame_request_callback_teardown(void)
{
    BObolViewController controller;
    BlockingFrameRequest state;
    controller.setFrameRequestCallback(blocking_frame_request, &state);

    std::thread dispatch([&controller]() {
	controller.markProgressiveWorkPending();
    });
    {
	std::unique_lock<std::mutex> lock(state.mutex);
	const bool entered = state.cv.wait_for(lock, std::chrono::seconds(2),
	    [&state]() { return state.entered; });
	if (!entered) {
	    state.release = true;
	    state.cv.notify_all();
	    lock.unlock();
	    dispatch.join();
	    CHECK(false, "frame callback dispatch starts");
	}
    }

    std::future<void> clear = std::async(std::launch::async, [&controller,
	&state]() {
	controller.clearFrameRequestCallback(&state);
    });
    const bool clear_waited = clear.wait_for(std::chrono::milliseconds(100)) ==
	std::future_status::timeout;
    {
	std::lock_guard<std::mutex> lock(state.mutex);
	state.release = true;
    }
    state.cv.notify_all();
    dispatch.join();
    const bool clear_finished = clear.wait_for(std::chrono::seconds(2)) ==
	std::future_status::ready;
    if (clear_finished)
	clear.get();

    CHECK(clear_waited && clear_finished && state.returned,
	  "callback teardown waits for an in-flight progressive dispatch");
    controller.clearProgressiveWorkPending();
    return 0;
}

static int
test_factory_endpoint_close_during_frame_request(void)
{
    FactoryTestState state;
    state.block_frames = true;
    struct bobol_host_factory factory = factory_test_desc(
	"endpoint-test-close-during-frame",
	10, BOBOL_HOST_CAP_PIXEL_PRESENT, &state);
    bobol_host_factory_token_t *token =
	bobol_host_factory_register(&factory);
    CHECK(token != NULL, "active-frame endpoint factory registers");

    struct bobol_host_desc desc = {};
    desc.struct_size = sizeof(desc);
    desc.mode = BOBOL_HOST_MODE_HEADLESS;
    desc.width = 16;
    desc.height = 12;
    desc.device_pixel_ratio = 1.0;
    desc.required_capabilities = BOBOL_HOST_CAP_PIXEL_PRESENT;

    bobol_display_endpoint_t *endpoint =
	bobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint != NULL, "active-frame endpoint creates");
    BObolViewController *controller =
	static_cast<BObolViewController *>(
	    bobol_display_endpoint_controller(endpoint));
    controller->clearRenderRequest();
    controller->clearProgressiveWorkPending();
    CHECK(bobol_display_endpoint_host_open(endpoint,
	  "endpoint-test-close-during-frame", &desc),
	  "active-frame endpoint opens its factory host");

    std::thread frame_request([controller]() {
	controller->markProgressiveWorkPending();
    });
    {
	std::unique_lock<std::mutex> lock(state.frame_mutex);
	const bool entered = state.frame_cv.wait_for(lock,
	    std::chrono::seconds(2), [&state]() {
		return state.frame_entered;
	    });
	if (!entered) {
	    state.frame_release = true;
	    state.frame_cv.notify_all();
	    lock.unlock();
	    frame_request.join();
	    bobol_display_endpoint_destroy(endpoint);
	    CHECK(false, "factory frame callback starts before endpoint teardown");
	}
    }

    std::future<void> close = std::async(std::launch::async, [endpoint]() {
	bobol_display_endpoint_destroy(endpoint);
    });
    const bool close_waited =
	close.wait_for(std::chrono::milliseconds(100)) ==
	std::future_status::timeout;
    {
	std::lock_guard<std::mutex> lock(state.frame_mutex);
	state.frame_release = true;
    }
    state.frame_cv.notify_all();
    frame_request.join();
    const bool close_finished = close.wait_for(std::chrono::seconds(2)) ==
	std::future_status::ready;
    if (close_finished)
	close.get();

    CHECK(close_waited && close_finished && state.frame_returned &&
	  state.closes == 1 && state.detaches == 1 && state.destroys == 1,
	  "endpoint close drains an active factory frame callback before host destroy");
    CHECK(bobol_host_factory_unregister(token),
	  "active-frame endpoint releases its factory after close");
    return 0;
}

static int
test_display_endpoint_contract(void)
{
    bobol_display_endpoint_t *endpoint =
	bobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint != NULL, "display endpoint creates an owned controller");
    CHECK(bobol_display_endpoint_controller(endpoint) != NULL,
	  "display endpoint exposes its controller");
    CHECK(bobol_display_endpoint_render_engine_get(endpoint) ==
	  BOBOL_RENDER_ENGINE_AUTO,
	  "display endpoint defaults to automatic renderer selection");
    CHECK(bobol_display_endpoint_render_engine_resolved_get(endpoint) ==
	  BOBOL_RENDER_ENGINE_AUTO,
	  "unbound automatic renderer remains explicitly unresolved");
    CHECK(bobol_display_endpoint_render_engine_set(endpoint,
	  BOBOL_RENDER_ENGINE_SW),
	  "display endpoint accepts a typed renderer policy");
    CHECK(bobol_display_endpoint_render_engine_get(endpoint) ==
	  BOBOL_RENDER_ENGINE_SW,
	  "display endpoint retains renderer policy");
    CHECK(bobol_display_endpoint_render_engine_capabilities(endpoint) ==
	  (BOBOL_RENDER_ENGINE_CAP_AUTO | BOBOL_RENDER_ENGINE_CAP_HW |
	   BOBOL_RENDER_ENGINE_CAP_SW | BOBOL_RENDER_ENGINE_CAP_RT |
	   BOBOL_RENDER_ENGINE_CAP_NONE |
	   BOBOL_RENDER_ENGINE_CAP_DIAGNOSTIC),
	  "unbound endpoint reports all selectable pre-host renderer policies");
    CHECK(!bobol_display_endpoint_render_engine_set(endpoint,
	  (enum bobol_render_engine)99),
	  "display endpoint rejects an invalid renderer policy");

    CHECK(bobol_display_endpoint_property_count() >= 34,
	  "display endpoint exposes the typed property registry");
    int found_renderer_property = 0;
    int found_resolved_renderer_property = 0;
    int found_title_property = 0;
    int found_visibility_property = 0;
    int found_vsync_property = 0;
    int found_perspective_property = 0;
    int found_faceplate_grid_property = 0;
    int found_framebuffer_composition_property = 0;
    int found_antialiasing_property = 0;
    int found_faceplate_font_properties = 0;
    int found_faceplate_color_properties = 0;
    for (size_t i = 0; i < bobol_display_endpoint_property_count(); i++) {
	struct bv_display_property_desc property = {};
	property.struct_size = sizeof(property);
	CHECK(bobol_display_endpoint_property_descriptor(i, &property) ==
	      BV_DISPLAY_PROPERTY_OK && property.name,
	      "display endpoint property descriptors are readable");
	if (bu_strcmp(property.name, "endpoint.renderer") == 0) {
	    found_renderer_property = property.type ==
		BV_DISPLAY_PROPERTY_ENUM &&
		(property.access & BV_DISPLAY_PROPERTY_WRITE) &&
		property.allowed_values &&
		bu_strcmp(property.allowed_values,
		    "auto,hw,sw,rt,none,diagnostic") == 0;
	}
	if (bu_strcmp(property.name, "endpoint.renderer.resolved") == 0) {
	    found_resolved_renderer_property = property.type ==
		BV_DISPLAY_PROPERTY_ENUM &&
		property.access == BV_DISPLAY_PROPERTY_READ &&
		property.allowed_values &&
		bu_strcmp(property.allowed_values,
		    "auto,hw,sw,rt,none,diagnostic") == 0;
	}
	if (bu_strcmp(property.name, "endpoint.title") == 0)
	    found_title_property = property.type ==
		BV_DISPLAY_PROPERTY_STRING &&
		property.required_host_capabilities == BOBOL_HOST_CAP_TOPLEVEL;
	if (bu_strcmp(property.name, "endpoint.visible") == 0)
	    found_visibility_property = property.type ==
		BV_DISPLAY_PROPERTY_BOOL &&
		property.required_host_capabilities == BOBOL_HOST_CAP_TOPLEVEL;
	if (bu_strcmp(property.name, "endpoint.vsync") == 0)
	    found_vsync_property = property.type ==
		BV_DISPLAY_PROPERTY_BOOL &&
		property.required_host_capabilities ==
		    BOBOL_HOST_CAP_PRESENT_VSYNC;
	if (bu_strcmp(property.name, "view.perspective") == 0)
	    found_perspective_property = property.type ==
		BV_DISPLAY_PROPERTY_DOUBLE &&
		(property.access & BV_DISPLAY_PROPERTY_WRITE) &&
		fabs(property.minimum) < 0.0001 && property.maximum < 180.0;
	if (bu_strcmp(property.name, "view.faceplate.grid.visible") == 0)
	    found_faceplate_grid_property = property.type ==
		BV_DISPLAY_PROPERTY_BOOL &&
		(property.access & BV_DISPLAY_PROPERTY_WRITE);
	if (bu_strcmp(property.name, "composition.framebuffer.mode") == 0)
	    found_framebuffer_composition_property = property.type ==
		BV_DISPLAY_PROPERTY_ENUM &&
		(property.access & BV_DISPLAY_PROPERTY_WRITE) &&
		property.allowed_values &&
		bu_strcmp(property.allowed_values, "off,overlay,underlay,interlay") == 0;
	if (bu_strcmp(property.name, "renderer.antialiasing") == 0)
	    found_antialiasing_property = property.type ==
		BV_DISPLAY_PROPERTY_BOOL &&
		(property.access & BV_DISPLAY_PROPERTY_WRITE);
	if (bu_strcmp(property.name, "view.faceplate.params.font_size") == 0 ||
	    bu_strcmp(property.name, "view.faceplate.center_dot.font_size") == 0 ||
	    bu_strcmp(property.name, "view.faceplate.scale.font_size") == 0) {
	    if (property.type == BV_DISPLAY_PROPERTY_UINT &&
		(property.access & BV_DISPLAY_PROPERTY_WRITE) &&
		fabs(property.minimum - 5.0) < 0.0001 &&
		fabs(property.maximum - 96.0) < 0.0001)
		found_faceplate_font_properties++;
	}
	if (bu_strcmp(property.name, "view.faceplate.params.color") == 0 ||
	    bu_strcmp(property.name, "view.faceplate.center_dot.color") == 0 ||
	    bu_strcmp(property.name, "view.faceplate.scale.color") == 0) {
	    if (property.type == BV_DISPLAY_PROPERTY_COLOR3 &&
		(property.access & BV_DISPLAY_PROPERTY_WRITE) &&
		fabs(property.minimum) < 0.0001 &&
		fabs(property.maximum - 1.0) < 0.0001)
		found_faceplate_color_properties++;
	}
    }
    CHECK(found_renderer_property,
	  "renderer property declares its type, access, and allowed values");
    CHECK(found_resolved_renderer_property,
	  "resolved renderer is a read-only typed endpoint property");
    CHECK(found_title_property && found_visibility_property,
	  "toplevel host properties declare typed capability predicates");
    CHECK(found_vsync_property,
	  "vsync property declares its presentation capability predicate");
    CHECK(found_perspective_property,
	  "perspective property declares a bounded double camera policy");
    CHECK(found_faceplate_grid_property,
	  "faceplate visibility is a writable endpoint policy");
    CHECK(found_framebuffer_composition_property,
	  "framebuffer composition declares its typed ordering values");
    CHECK(found_antialiasing_property,
	  "antialiasing is a writable renderer policy");
    CHECK(found_faceplate_font_properties == 3,
	  "faceplate font policies declare bounded writable uint values");
    CHECK(found_faceplate_color_properties == 3,
	  "faceplate color policies declare bounded writable color values");

    struct bv_display_property_value property_value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    property_value.type = BV_DISPLAY_PROPERTY_STRING;
    property_value.string_value = "No host";
    CHECK(bobol_display_endpoint_property_set(endpoint, "endpoint.title",
	  &property_value) == BV_DISPLAY_PROPERTY_UNSUPPORTED,
	  "host properties fail explicitly without a compatible toplevel host");
    property_value_init(&property_value);
    CHECK(bobol_display_endpoint_property_get(endpoint,
	  "endpoint.renderer", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK &&
	  property_value.type == BV_DISPLAY_PROPERTY_ENUM &&
	  bu_strcmp(property_value.string_value, "sw") == 0,
	  "typed renderer property reads endpoint policy");
    property_value_init(&property_value);
    CHECK(bobol_display_endpoint_property_get(endpoint,
	  "endpoint.renderer.supported", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK && property_value.string_value &&
	  bu_strcmp(property_value.string_value,
	      "auto,hw,sw,rt,none,diagnostic") == 0,
	  "typed renderer support property reports the current binding");
    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_ENUM;
    property_value.string_value = "quality";
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "controller.software_wire", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK,
	  "typed software wire property accepts QUALITY");
    BObolViewController *endpoint_controller =
	static_cast<BObolViewController *>(
	    bobol_display_endpoint_controller(endpoint));
    CHECK(endpoint_controller->getSceneRoot() &&
	  endpoint_controller->getSceneRoot()->isOfType(
	      SoBRLSceneGroup::getClassTypeId()),
	  "endpoint-created controller starts with an authoritative scene root");
    CHECK(endpoint_controller->getSoftwareWireMode() ==
	  BObolViewController::SOFTWARE_WIRE_QUALITY,
	  "typed software wire property updates the controller");

    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_BOOL;
    property_value.bool_value = 0;
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "renderer.depth_test", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK &&
	  !endpoint_controller->isDepthTestEnabled(),
	  "typed depth property updates the Obol render environment");
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "renderer.lighting", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK &&
	  !endpoint_controller->isLightingEnabled(),
	  "typed lighting property updates the Obol render environment");
    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_ENUM;
    property_value.string_value = "mged";
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "renderer.lighting.profile", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK &&
	  endpoint_controller->getLightingProfile() ==
	      BObolViewController::LIGHTING_MGED,
	  "typed lighting profile selects the historical MGED rig");
    property_value_init(&property_value);
    CHECK(bobol_display_endpoint_property_get(endpoint,
	  "renderer.lighting.profile", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK && property_value.string_value &&
	  bu_strcmp(property_value.string_value, "mged") == 0,
	  "typed lighting profile round trips render policy");
    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_ENUM;
    property_value.string_value = "studio";
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "renderer.lighting.profile", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK &&
	  endpoint_controller->getLightingProfile() ==
	      BObolViewController::LIGHTING_STUDIO,
	  "typed lighting profile restores the studio rig");
    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_COLOR3;
    property_value.color3[0] = 0.25;
    property_value.color3[1] = 0.50;
    property_value.color3[2] = 0.75;
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "renderer.headlight.color", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK &&
	  endpoint_controller->getHeadlightColor() ==
	      SbColor(0.25f, 0.50f, 0.75f),
	  "typed headlight color updates the common render environment");
    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_DOUBLE;
    property_value.double_value = 0.4;
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "renderer.headlight.intensity", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK &&
	  std::fabs(endpoint_controller->getHeadlightIntensity() - 0.4f) <
	      1.0e-6f,
	  "typed headlight intensity updates the common render environment");
    property_value_init(&property_value);
    CHECK(bobol_display_endpoint_property_get(endpoint,
	  "renderer.headlight.color", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK &&
	  std::fabs(property_value.color3[0] - 0.25) < 1.0e-9 &&
	  std::fabs(property_value.color3[1] - 0.50) < 1.0e-9 &&
	  std::fabs(property_value.color3[2] - 0.75) < 1.0e-9,
	  "typed headlight color round trips render policy");
    endpoint_controller->setHeadlightColor(SbColor(1.0f, 1.0f, 1.0f));
    endpoint_controller->setHeadlightIntensity(1.0f);
    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_BOOL;
    property_value.bool_value = 0;
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "renderer.transparency", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK &&
	  !endpoint_controller->isTransparencyEnabled(),
	  "typed transparency property updates the Obol render action");
    property_value_init(&property_value);
    CHECK(bobol_display_endpoint_property_get(endpoint,
	  "renderer.transparency", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK && !property_value.bool_value,
	  "typed transparency property round trips render policy");
    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_BOOL;
    property_value.bool_value = 1;
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "renderer.antialiasing", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK &&
	  endpoint_controller->isAntialiasingEnabled(),
	  "typed antialiasing property updates the Obol render action");
    property_value_init(&property_value);
    CHECK(bobol_display_endpoint_property_get(endpoint,
	  "renderer.antialiasing", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK && property_value.bool_value,
	  "typed antialiasing property round trips render policy");
    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_BOOL;
    property_value.bool_value = 2;
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "renderer.antialiasing", &property_value) ==
	  BV_DISPLAY_PROPERTY_INVALID,
	  "typed Boolean policy rejects values outside its declared range");
    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_DOUBLE;
    property_value.double_value = -4.0;
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "renderer.clip.minimum", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK,
	  "typed clip minimum updates controller camera policy");
    property_value.double_value = 5.0;
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "renderer.clip.maximum", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK,
	  "typed clip maximum updates controller camera policy");
    double clip_minimum = 0.0;
    double clip_maximum = 0.0;
    endpoint_controller->getClipBounds(clip_minimum, clip_maximum);
    CHECK(fabs(clip_minimum + 4.0) < 0.0001 &&
	  fabs(clip_maximum - 5.0) < 0.0001,
	  "typed clip bounds preserve the ordered controller range");
    struct bv_context *clip_view = bv_context_create();
    CHECK(clip_view != NULL &&
	  bv_context_dimensions_set(clip_view, 320, 240) &&
	  bv_zclip_set(bv_context_view(clip_view), 1) &&
	  endpoint_controller->syncCameraFromViewContext(clip_view),
	  "controller synchronizes enabled retained clipping");
    SoClipPlane *clip_minimum_node = find_clip_plane(
	endpoint_controller->getRenderRoot(), "BObolClipMinimum");
    SoClipPlane *clip_maximum_node = find_clip_plane(
	endpoint_controller->getRenderRoot(), "BObolClipMaximum");
    CHECK(clip_minimum_node && clip_maximum_node &&
	  clip_minimum_node->on.getValue() &&
	  clip_maximum_node->on.getValue(),
	  "zclip enables the retained Obol clipping-plane pair");
    CHECK(bv_zclip_set(bv_context_view(clip_view), 0) &&
	  endpoint_controller->syncCameraFromViewContext(clip_view) &&
	  !clip_minimum_node->on.getValue() &&
	  !clip_maximum_node->on.getValue(),
	  "disabling zclip preserves but deactivates retained clip planes");
    bv_context_destroy(clip_view);
    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_BOOL;
    property_value.bool_value = 1;
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "renderer.depth_cue", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK &&
	  endpoint_controller->isDepthCueEnabled(),
	  "typed depth-cue property updates the Obol render environment");

    CHECK(bobol_display_endpoint_property_get(endpoint, "view.zclip",
	  &property_value) == BV_DISPLAY_PROPERTY_UNSUPPORTED,
	  "external view property reports unsupported without an owner");
    CHECK(bobol_display_endpoint_property_get(endpoint, "view.perspective",
	  &property_value) == BV_DISPLAY_PROPERTY_UNSUPPORTED,
	  "external perspective property reports unsupported without an owner");
    CHECK(bobol_display_endpoint_property_get(endpoint,
	  "composition.framebuffer.mode", &property_value) ==
	  BV_DISPLAY_PROPERTY_UNSUPPORTED,
	  "external composition property reports unsupported without an owner");
    TestPropertyProvider provider = {};
    CHECK(bobol_display_endpoint_property_provider_set(endpoint,
	  test_property_provider_get, test_property_provider_set, &provider),
	  "endpoint accepts an external property owner");
    int faceplate_visible = -1;
    CHECK(bobol_display_endpoint_input_faceplate_toggle_apply(endpoint,
	  NULL, BOBOL_ACTION_TOGGLE_ADC, &faceplate_visible) &&
	  faceplate_visible == 1 && provider.adcVisible == 1,
	  "shared faceplate action uses the endpoint ADC property owner");
    CHECK(bobol_display_endpoint_input_faceplate_toggle_apply(endpoint,
	  NULL, BOBOL_ACTION_TOGGLE_MODEL_AXES, &faceplate_visible) &&
	  faceplate_visible == 1 && provider.modelAxesVisible == 1,
	  "shared faceplate action uses the endpoint model-axes property owner");
    CHECK(bobol_display_endpoint_input_faceplate_toggle_apply(endpoint,
	  NULL, BOBOL_ACTION_TOGGLE_VIEW_AXES, &faceplate_visible) &&
	  faceplate_visible == 1 && provider.viewAxesVisible == 1,
	  "shared faceplate action uses the endpoint view-axes property owner");
    property_value_init(&property_value);
    property_value.bool_value = 1;
    CHECK(bobol_display_endpoint_property_set(endpoint, "view.zclip",
	  &property_value) == BV_DISPLAY_PROPERTY_OK && provider.zclip == 1,
	  "external property setter updates owner state");
    property_value.bool_value = 0;
    CHECK(bobol_display_endpoint_property_get(endpoint, "view.zclip",
	  &property_value) == BV_DISPLAY_PROPERTY_OK &&
	  property_value.bool_value == 1,
	  "external property getter reads owner state");
    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_DOUBLE;
    property_value.double_value = 45.0;
    CHECK(bobol_display_endpoint_property_set(endpoint, "view.perspective",
	  &property_value) == BV_DISPLAY_PROPERTY_OK &&
	  fabs(provider.perspective - 45.0) < 0.0001,
	  "external perspective setter updates owner state");
    property_value_init(&property_value);
    CHECK(bobol_display_endpoint_property_get(endpoint, "view.perspective",
	  &property_value) == BV_DISPLAY_PROPERTY_OK &&
	  property_value.type == BV_DISPLAY_PROPERTY_DOUBLE &&
	  fabs(property_value.double_value - 45.0) < 0.0001,
	  "external perspective getter reads owner state");
    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_ENUM;
    property_value.string_value = "interlay";
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "composition.framebuffer.mode", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK && provider.framebufferMode == 3,
	  "external composition setter updates owner state");
    property_value_init(&property_value);
    CHECK(bobol_display_endpoint_property_get(endpoint,
	  "composition.framebuffer.mode", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK &&
	  property_value.type == BV_DISPLAY_PROPERTY_ENUM &&
	  property_value.string_value &&
	  bu_strcmp(property_value.string_value, "interlay") == 0,
	  "external composition getter reads owner state");
    CHECK(bobol_display_endpoint_property_provider_set(endpoint,
	  NULL, NULL, NULL),
	  "endpoint clears its external property owner");

    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_COLOR3;
    property_value.color3[0] = 0.125;
    property_value.color3[1] = 0.25;
    property_value.color3[2] = 0.5;
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "controller.background.bottom", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK,
	  "typed background property accepts normalized RGB");
    property_value_init(&property_value);
    CHECK(bobol_display_endpoint_property_get(endpoint,
	  "controller.background.bottom", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK &&
	  fabs(property_value.color3[0] - 0.125) < 0.0001 &&
	  fabs(property_value.color3[1] - 0.25) < 0.0001 &&
	  fabs(property_value.color3[2] - 0.5) < 0.0001,
	  "typed background property round trips controller color");

    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_UINT;
    property_value.uint_value = 320;
    CHECK(bobol_display_endpoint_property_set(endpoint, "endpoint.width",
	  &property_value) == BV_DISPLAY_PROPERTY_OK,
	  "typed width property resizes the endpoint");
    property_value.uint_value = 240;
    CHECK(bobol_display_endpoint_property_set(endpoint, "endpoint.height",
	  &property_value) == BV_DISPLAY_PROPERTY_OK,
	  "typed height property resizes the endpoint");
    const SbVec2s property_viewport =
	endpoint_controller->getViewportRegion().getWindowSize();
    CHECK(property_viewport[0] == 320 && property_viewport[1] == 240,
	  "typed size properties update the controller viewport");

    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_STRING;
    property_value.string_value = "invalid";
    CHECK(bobol_display_endpoint_property_set(endpoint, "endpoint.host",
	  &property_value) == BV_DISPLAY_PROPERTY_READ_ONLY,
	  "typed property API rejects writes to read-only host identity");
    CHECK(bobol_display_endpoint_property_get(endpoint,
	  "does.not.exist", &property_value) ==
	  BV_DISPLAY_PROPERTY_UNKNOWN,
	  "typed property API distinguishes unknown properties");
    property_value.type = BV_DISPLAY_PROPERTY_DOUBLE;
    CHECK(bobol_display_endpoint_property_set(endpoint, "endpoint.width",
	  &property_value) == BV_DISPLAY_PROPERTY_INVALID,
	  "typed property API rejects mismatched value types");

    property_value_init(&property_value);
    property_value.type = BV_DISPLAY_PROPERTY_ENUM;
    property_value.string_value = "rt";
    CHECK(bobol_display_endpoint_property_set(endpoint,
	  "endpoint.renderer", &property_value) ==
	  BV_DISPLAY_PROPERTY_OK &&
	  bobol_display_endpoint_render_engine_get(endpoint) ==
	  BOBOL_RENDER_ENGINE_RT,
	  "typed renderer property selects retained librt rendering");
    CHECK(bobol_display_endpoint_render_engine_set(endpoint,
	  BOBOL_RENDER_ENGINE_SW),
	  "software renderer policy remains selectable after retained librt");

    BObolWindowHost borrowed_host;
    CHECK(!bobol_display_endpoint_host_bind(endpoint, &borrowed_host, 0),
	  "explicit software renderer rejects an untyped direct host");
    CHECK(bobol_display_endpoint_render_engine_set(endpoint,
	  BOBOL_RENDER_ENGINE_AUTO),
	  "automatic renderer policy permits custom direct host binding");
    CHECK(bobol_display_endpoint_host_bind(endpoint, &borrowed_host, 0),
	  "display endpoint binds a borrowed host");
    CHECK(bobol_display_endpoint_render_engine_resolved_get(endpoint) ==
	  BOBOL_RENDER_ENGINE_AUTO,
	  "custom direct host does not invent renderer capabilities");
    CHECK(bobol_display_endpoint_host(endpoint) == &borrowed_host,
	  "display endpoint reports its host");
    CHECK(bobol_display_endpoint_framebuffer_window_host(endpoint) ==
	  &borrowed_host,
	  "direct endpoint host remains available for framebuffer attachment");
    CHECK(borrowed_host.getController() ==
	  bobol_display_endpoint_controller(endpoint),
	  "bound host borrows the endpoint controller");
    bobol_display_endpoint_host_detach(endpoint);
    CHECK(bobol_display_endpoint_host(endpoint) == NULL &&
	  borrowed_host.getController() == NULL,
	  "display endpoint detaches a borrowed host");

    CountingWindowHost::destroyed = 0;
    CountingWindowHost *owned_host = new CountingWindowHost;
    CHECK(bobol_display_endpoint_host_bind(endpoint, owned_host,
	  BOBOL_ENDPOINT_OWN_HOST),
	  "display endpoint binds an owned host");
    CHECK(bobol_display_endpoint_host_bind(endpoint, owned_host, 0),
	  "display endpoint permits idempotent host binding");
    bobol_display_endpoint_destroy(endpoint);
    CHECK(CountingWindowHost::destroyed == 1,
	  "idempotent binding does not relinquish owned host lifetime");

    BObolViewController borrowed_controller;
    BObolWindowHost wrapper_host;
    {
	BObolDisplayEndpoint wrapper(&borrowed_controller, false);
	CHECK(wrapper.isValid() && wrapper.controller() == &borrowed_controller,
	      "C++ endpoint wrapper borrows an explicit controller");
	CHECK(wrapper.bindHost(&wrapper_host) && wrapper.host() == &wrapper_host,
	      "C++ endpoint wrapper binds a typed host");
	CHECK(!wrapper.setRenderEngine(BOBOL_RENDER_ENGINE_HW),
	      "explicit hardware policy rejects an untyped direct host");
	CHECK(!wrapper.supportsRenderEngine(BOBOL_RENDER_ENGINE_HW) &&
	      wrapper.supportsRenderEngine(BOBOL_RENDER_ENGINE_NONE) &&
	      (wrapper.renderEngineCapabilities() &
	       BOBOL_RENDER_ENGINE_CAP_DIAGNOSTIC),
	      "C++ endpoint wrapper exposes binding-specific renderer capabilities");
	CHECK(wrapper.setRenderEngine(BOBOL_RENDER_ENGINE_NONE) &&
	      wrapper.renderEngine() == BOBOL_RENDER_ENGINE_NONE &&
	      wrapper.resolvedRenderEngine() == BOBOL_RENDER_ENGINE_NONE,
	      "C++ endpoint wrapper exposes renderer policy");
    }
    CHECK(wrapper_host.getController() == NULL &&
	  borrowed_controller.getRenderRoot() != NULL,
	  "endpoint destruction detaches its host and restores a borrowed controller after none");
    return 0;
}

static int
test_display_endpoint_input(void)
{
    bobol_display_endpoint_t *first =
	bobol_display_endpoint_create(NULL, 0);
    bobol_display_endpoint_t *second =
	bobol_display_endpoint_create(NULL, 0);
    CHECK(first && second, "display endpoints allocate input contexts");

    InputActionState state;
    BObolInputEvent event;
    event.type = BOBOL_INPUT_KEY_PRESS;
    event.key = 'F';
    CHECK(bobol_display_endpoint_input_action_handler_set(first,
	  test_input_action, &state) &&
	  bobol_display_endpoint_input_dispatch(first, &event) == 1 &&
	  state.action == BOBOL_ACTION_VIEW_FRONT,
	  "display endpoint dispatches its default input profile");
    CHECK(bobol_display_endpoint_input_dispatch(second, &event) == 0,
	  "display endpoint with no action handler leaves input unhandled");

    BObolInputBinding semantic_binding;
    semantic_binding.eventType = BOBOL_INPUT_POINTER_RELEASE;
    semantic_binding.key = BOBOL_INPUT_ANY;
    semantic_binding.button = 0;
    semantic_binding.requiredModifiers = BOBOL_INPUT_MOD_NONE;
    semantic_binding.forbiddenModifiers = BOBOL_INPUT_MOD_SHIFT |
	BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_ALT |
	BOBOL_INPUT_MOD_META;
    semantic_binding.priority = 0;
    semantic_binding.action = TEST_LAYER_SELECT_COMMIT;
    InputActionState semantic_state;
    semantic_state.result = BOBOL_INPUT_RESULT_DEFERRED;
    const BObolInputActionLayer semantic_layer = {"endpoint-selection",
	&semantic_binding, 1, test_input_action};
    CHECK(bobol_display_endpoint_input_action_layer_set(first,
	  &semantic_layer, &semantic_state, &semantic_state),
	  "display endpoint installs a scoped application action layer");
    event.type = BOBOL_INPUT_POINTER_RELEASE;
    event.key = BOBOL_INPUT_ANY;
    event.button = 0;
    event.modifiers = BOBOL_INPUT_MOD_NONE;
    const int view_calls_before_layer = state.calls;
    CHECK(bobol_display_endpoint_input_dispatch(first, &event) ==
	  BOBOL_INPUT_RESULT_DEFERRED && semantic_state.calls == 1 &&
	  semantic_state.action == TEST_LAYER_SELECT_COMMIT &&
	  state.calls == view_calls_before_layer,
	  "display endpoint preserves layer results without replacing view handling");
    semantic_state.result = BOBOL_INPUT_RESULT_UNHANDLED;
    CHECK(bobol_display_endpoint_input_dispatch(first, &event) ==
	  BOBOL_INPUT_RESULT_HANDLED && semantic_state.calls == 2 &&
	  state.calls == view_calls_before_layer + 1 &&
	  state.action == BOBOL_ACTION_VIEW_PRIMARY_RELEASE,
	  "display endpoint falls through an unhandled application action");
    CHECK(bobol_display_endpoint_input_action_layer_clear_if(first,
	  &semantic_state),
	  "display endpoint conditionally clears its application layer");

    BObolInputBinding binding;
    binding.eventType = BOBOL_INPUT_KEY_PRESS;
    binding.key = 'P';
    binding.button = BOBOL_INPUT_ANY;
    binding.requiredModifiers = BOBOL_INPUT_MOD_NONE;
    binding.forbiddenModifiers = BOBOL_INPUT_MOD_NONE;
    binding.priority = 0;
    binding.action = BOBOL_ACTION_VIEW_PAN;
    BObolInputProfile profile = {"endpoint-custom", &binding, 1};
    CHECK(bobol_display_endpoint_input_profile_set(first, &profile) &&
	  bobol_display_endpoint_input_dispatch(first, &event) == 0,
	  "endpoint-local profiles replace only their own bindings");
    event.type = BOBOL_INPUT_KEY_PRESS;
    event.key = 'P';
    CHECK(bobol_display_endpoint_input_dispatch(first, &event) == 1 &&
	  state.action == BOBOL_ACTION_VIEW_PAN,
	  "endpoint dispatches a custom semantic profile");

    bobol_display_endpoint_destroy(second);
    bobol_display_endpoint_destroy(first);
    return 0;
}

static int
texture_matches(SoTexture2 *texture, int width, int height, int channels)
{
    if (!texture)
	return 0;
    int texWidth = 0;
    int texHeight = 0;
    int texChannels = 0;
    const unsigned char *pixels = texture->getImageData(texWidth, texHeight, texChannels);
    return pixels && texWidth == width && texHeight == height && texChannels == channels;
}

static int
test_imgstream_display_host_bridge(void)
{
    BObolWindowHost host;

    CHECK(bobol_window_host_open_display_framebuffer(NULL, "/dev/qtgl", 5, 4) == NULL,
	  "null Obol display host is rejected");

    imgstream_fb_t *fb = bobol_window_host_open_display_framebuffer(&host, "/dev/qtgl", 5, 4);
    CHECK(fb != NULL, "display framebuffer opens through Obol host");
    CHECK(host.isOpen(), "display framebuffer opens host");
    CHECK(host.getFramebufferCount() == 1, "host records one framebuffer");

    SoBRLImageSource *source = host.getFramebufferImageSource(fb);
    SoBRLViewportImage *viewport = host.getFramebufferViewportImage(fb);
    CHECK(source != NULL, "display host creates image source");
    CHECK(viewport != NULL, "display host creates viewport image");
    CHECK(source->getStream() == imgstream_fb_stream(fb),
	  "display host image source borrows framebuffer stream");
    CHECK(viewport->getImageSource() == source,
	  "display host viewport image references source");
    CHECK(viewport->getTextureNode() != NULL &&
	  texture_matches(viewport->getTextureNode(), 5, 4, 3),
	  "display host realizes framebuffer texture");
    SoGroup *underlay = host.getController()->getFramebufferUnderlayRoot();
    SoGroup *interlay = host.getController()->getFramebufferInterlayRoot();
    SoGroup *overlay = host.getController()->getFramebufferOverlayRoot();

    CHECK(underlay && interlay && overlay && overlay->findChild(viewport) >= 0 &&
	  underlay->findChild(viewport) < 0 &&
	  interlay->findChild(viewport) < 0 &&
	  viewport->layer.getValue() == SoBRLViewportImage::OVERLAY,
	  "display host attaches a framebuffer image to the retained overlay layer");
    CHECK(host.setFramebufferComposition(fb,
	  BOBOL_FRAMEBUFFER_COMPOSITION_UNDERLAY) == 0 &&
	  underlay->findChild(viewport) >= 0 && overlay->findChild(viewport) < 0 &&
	  viewport->visible.getValue() == TRUE &&
	  viewport->layer.getValue() == SoBRLViewportImage::UNDERLAY,
	  "display host moves framebuffer image beneath the CAD render batch");
    CHECK(host.setFramebufferComposition(fb,
	  BOBOL_FRAMEBUFFER_COMPOSITION_INTERLAY) == 0 &&
	  interlay->findChild(viewport) >= 0 &&
	  underlay->findChild(viewport) < 0 && overlay->findChild(viewport) < 0 &&
	  viewport->visible.getValue() == TRUE &&
	  viewport->layer.getValue() == SoBRLViewportImage::INTERLAY,
	  "display host moves framebuffer image between CAD and view-local features");
    CHECK(host.setFramebufferComposition(fb,
	  BOBOL_FRAMEBUFFER_COMPOSITION_OFF) == 0 &&
	  underlay->findChild(viewport) < 0 &&
	  interlay->findChild(viewport) < 0 && overlay->findChild(viewport) < 0 &&
	  viewport->visible.getValue() == FALSE,
	  "display host removes disabled framebuffer image from retained composition");
    CHECK(host.setFramebufferComposition(fb,
	  BOBOL_FRAMEBUFFER_COMPOSITION_OVERLAY) == 0 &&
	  overlay->findChild(viewport) >= 0 && viewport->visible.getValue() == TRUE,
	  "display host restores framebuffer image above the CAD render batch");
    (void)host.getController()->consumeRenderRequest(NULL);
    CHECK(host.setFramebufferComposition(fb,
	  BOBOL_FRAMEBUFFER_COMPOSITION_OVERLAY) == 0 &&
	  !host.getController()->isRenderRequested() &&
	  overlay->findChild(viewport) >= 0,
	  "repeating a framebuffer composition mode is an idle no-op");

    unsigned char red[3] = {255, 0, 0};
    CHECK(imgstream_fb_writerect(fb, 2, 1, 1, 1, red) == 1,
	  "framebuffer write updates stream");
    CHECK(imgstream_fb_flush(fb) == 0, "display host flush syncs stream");
    CHECK(viewport->realizedDirtyRevision.getValue() == source->dirtyRevision.getValue(),
	  "display host flush refreshes viewport dirty revision");

    CHECK(imgstream_fb_view(fb, 3, 2, 4, 4) == 0,
	  "display host accepts framebuffer view transform");
    CHECK(float_equal(viewport->sourceCenter.getValue()[0], 3.0f) &&
	  float_equal(viewport->sourceCenter.getValue()[1], 2.0f) &&
	  float_equal(viewport->sourceZoom.getValue(), 4.0f),
	  "display host maps framebuffer view to viewport image");

    CHECK(imgstream_fb_cursor(fb, 1, 4, 3) == 0,
	  "display host accepts cursor state");
    CHECK(viewport->cursorVisible.getValue() == TRUE &&
	  float_equal(viewport->cursorImagePosition.getValue()[0], 4.0f) &&
	  float_equal(viewport->cursorImagePosition.getValue()[1], 3.0f),
	  "display host maps cursor state to viewport image");

    CHECK(imgstream_fb_viewport(fb, 1, 2, 11, 8) == 0,
	  "display host accepts viewport state");
    CHECK(host.getController()->getViewportRegion().getWindowSize()[0] == 10 &&
	  host.getController()->getViewportRegion().getWindowSize()[1] == 6,
	  "display host maps viewport state to controller size");

    unsigned char cursorBits[1] = {0xff};
    CHECK(imgstream_fb_setcursor(fb, cursorBits, 8, 1, 0, 0) == 0,
	  "display host accepts custom cursor shape");
    CHECK(viewport->cursorShape.getValue() == SoBRLViewportImage::CURSOR_CUSTOM,
	  "display host records custom cursor policy");

    CHECK(imgstream_fb_poll(fb) == 0, "display host poll succeeds");
    CHECK(imgstream_fb_poll_rate(fb) == 0, "display host poll rate defaults to zero");

    imgstream_fb_close(fb);
    CHECK(host.getFramebufferCount() == 0, "display host closes framebuffer attachment");

    CHECK(imgstream_fb_open("/dev/qtgl", 2, 2) == NULL,
	  "display framebuffer requires an explicit host");
    return 0;
}

static int
test_framebuffer_stream_helper(void)
{
    class RejectingFramebufferHost : public BObolWindowHost {
    public:
	int openFramebuffer(imgstream_fb_t *UNUSED(fb),
	    const imgstream_fb_spec_info_t *UNUSED(info)) override
	{
	    return -1;
	}
    };
    class ConditionalCompositionHost : public BObolWindowHost {
    public:
	ConditionalCompositionHost() : reject(false) { }

	int setFramebufferComposition(imgstream_fb_t *fb,
	    BObolFramebufferComposition composition) override
	{
	    return reject ? -1 :
		BObolWindowHost::setFramebufferComposition(fb, composition);
	}

	bool reject;
    };

    BObolWindowHost host;
    BObolWindowHost second_host;
    RejectingFramebufferHost rejecting_host;
    ConditionalCompositionHost composition_host;
    BObolFramebufferStream stream(&host);

    CHECK(stream.configure(4, 3) == 0, "framebuffer stream accepts size");

    BObolFramebufferInfo info;
    CHECK(stream.info(&info) == 0, "framebuffer stream reports dimensions");
    CHECK(info.width == 4 && info.height == 3 &&
	  info.max_width == 4 && info.max_height == 3,
	  "framebuffer stream dimensions match configuration");
    CHECK(stream.framebuffer() != NULL, "framebuffer stream opens imgstream framebuffer");
    CHECK(host.getFramebufferCount() == 1,
	  "framebuffer stream attaches one image to Obol host");

    imgstream_fb_t *fb = stream.framebuffer();
    SoBRLImageSource *source = host.getFramebufferImageSource(fb);
    SoBRLViewportImage *viewport = host.getFramebufferViewportImage(fb);
    CHECK(source != NULL && viewport != NULL,
	  "framebuffer stream exposes host image nodes");

    unsigned char blue[3] = {0, 0, 255};
    CHECK(stream.writerect(1, 1, 1, 1, blue) == 1,
	  "framebuffer stream writes pixels");
    CHECK(source->hasPendingStreamUpdate() == TRUE,
	  "framebuffer stream write marks source update pending");
    CHECK(stream.present() == 0,
	  "framebuffer stream presents pending pixels to Obol image");
    CHECK(viewport->realizedDirtyRevision.getValue() ==
	  source->dirtyRevision.getValue(),
	  "framebuffer stream present refreshes Obol image");
    CHECK(source->hasPendingStreamUpdate() == FALSE,
	  "framebuffer stream present clears pending source update");

    CHECK(stream.view(2, 1, 3, 3) == 0,
	  "framebuffer stream records view state");
    CHECK(stream.present() == 0,
	  "framebuffer stream presents view state");
    CHECK(float_equal(viewport->sourceCenter.getValue()[0], 2.0f) &&
	  float_equal(viewport->sourceCenter.getValue()[1], 1.0f) &&
	  float_equal(viewport->sourceZoom.getValue(), 3.0f),
	  "framebuffer stream maps view state to viewport image");

    CHECK(stream.cursor(1, 3, 2) == 0,
	  "framebuffer stream records cursor state");
    CHECK(stream.present() == 0,
	  "framebuffer stream presents cursor state");
    CHECK(viewport->cursorVisible.getValue() == TRUE &&
	  float_equal(viewport->cursorImagePosition.getValue()[0], 3.0f) &&
	  float_equal(viewport->cursorImagePosition.getValue()[1], 2.0f),
	  "framebuffer stream maps cursor state to viewport image");

    CHECK(stream.setComposition(BOBOL_FRAMEBUFFER_COMPOSITION_UNDERLAY) == 0 &&
	  stream.composition() == BOBOL_FRAMEBUFFER_COMPOSITION_UNDERLAY &&
	  host.getController()->getFramebufferUnderlayRoot()->findChild(viewport) >= 0 &&
	  host.getController()->getFramebufferOverlayRoot()->findChild(viewport) < 0,
	  "framebuffer stream owns and applies typed composition state");

    composition_host.reject = true;
    CHECK(stream.setHost(&composition_host,
	      BOBOL_FRAMEBUFFER_COMPOSITION_INTERLAY) == -1 &&
	  stream.host() == &host &&
	  stream.composition() == BOBOL_FRAMEBUFFER_COMPOSITION_UNDERLAY &&
	  stream.framebuffer() == fb && host.getFramebufferCount() == 1 &&
	  composition_host.getFramebufferCount() == 0 &&
	  host.getController()->getFramebufferUnderlayRoot()->findChild(viewport) >= 0,
	  "composition rejection leaves the old framebuffer host and policy live");

    CHECK(stream.setHost(&rejecting_host) == -1 && stream.host() == &host &&
	  stream.framebuffer() == fb && host.getFramebufferCount() == 1 &&
	  rejecting_host.getFramebufferCount() == 0,
	  "failed framebuffer rehost leaves the live attachment unchanged");

    CHECK(stream.setHost(&second_host) == 0 &&
	  stream.host() == &second_host && stream.framebuffer() == fb &&
	  host.getFramebufferCount() == 0 &&
	  second_host.getFramebufferCount() == 1,
	  "live framebuffer stream moves to a second host without reopening");
    SoBRLImageSource *second_source =
	second_host.getFramebufferImageSource(fb);
    SoBRLViewportImage *second_viewport =
	second_host.getFramebufferViewportImage(fb);
    unsigned char moved_pixel[3] = {0, 0, 0};
    CHECK(second_source != NULL && second_source != source &&
	  second_viewport != NULL && second_viewport != viewport &&
	  stream.read(1, 1, moved_pixel, 1) == 1 &&
	  moved_pixel[0] == blue[0] && moved_pixel[1] == blue[1] &&
	  moved_pixel[2] == blue[2],
	  "framebuffer rehost preserves pixels and replaces only retained host nodes");
    CHECK(stream.composition() == BOBOL_FRAMEBUFFER_COMPOSITION_UNDERLAY &&
	  second_host.getController()->getFramebufferUnderlayRoot()->findChild(
	      second_viewport) >= 0 &&
	  second_host.getController()->getFramebufferOverlayRoot()->findChild(
	      second_viewport) < 0,
	  "framebuffer rehost preserves its typed composition layer");
    CHECK(stream.present() == 0 &&
	  float_equal(second_viewport->sourceCenter.getValue()[0], 2.0f) &&
	  float_equal(second_viewport->sourceCenter.getValue()[1], 1.0f) &&
	  float_equal(second_viewport->sourceZoom.getValue(), 3.0f) &&
	  second_viewport->cursorVisible.getValue() == TRUE &&
	  float_equal(second_viewport->cursorImagePosition.getValue()[0], 3.0f) &&
	  float_equal(second_viewport->cursorImagePosition.getValue()[1], 2.0f),
	  "moved framebuffer republishes preserved view and cursor state");
    CHECK(stream.setHost(&host) == 0 && stream.framebuffer() == fb &&
	  host.getFramebufferCount() == 1 &&
	  second_host.getFramebufferCount() == 0,
	  "framebuffer stream moves back without leaving stale retained nodes");

    stream.close();
    CHECK(host.getFramebufferCount() == 0,
	  "framebuffer stream close detaches host image nodes");

    CHECK(stream.configure(2, 2) == 0 && stream.ensure() == 0,
	  "framebuffer stream reopens after close");
    CHECK(host.getFramebufferCount() == 1,
	  "framebuffer stream reattaches after close");
    imgstream_fb_t *detached_fb = stream.framebuffer();
    CHECK(stream.setHost(NULL) == 0 && stream.host() == NULL &&
	  stream.framebuffer() == detached_fb && host.getFramebufferCount() == 0,
	  "framebuffer stream detaches presentation without closing storage");
    CHECK(stream.setHost(&second_host) == 0 &&
	  stream.framebuffer() == detached_fb &&
	  second_host.getFramebufferCount() == 1,
	  "detached framebuffer stream reattaches with the same identity");
    stream.close();
    CHECK(second_host.getFramebufferCount() == 0,
	  "reattached framebuffer closes without stale host nodes");
    return 0;
}

static int
test_framebuffer_host_teardown(void)
{
    for (int i = 0; i < 16; i++) {
	BObolWindowHost *host = new BObolWindowHost;
	BObolFramebufferStream *stream = new BObolFramebufferStream(host);
	CHECK(stream->configure(3, 2) == 0 && stream->ensure() == 0,
	      "stream opens before host teardown");

	imgstream_fb_t *display = bobol_window_host_open_display_framebuffer(
	    host, "/dev/qtgl", 3, 2);
	CHECK(display != NULL, "display framebuffer opens before host teardown");
	unsigned char pixel[3] = {0, 127, 255};
	CHECK(stream->write(0, 0, pixel, 1) == 1 && stream->present() == 0,
	      "stream presents an active image before host teardown");
	imgstream_fb_t *stream_fb = stream->framebuffer();
	CHECK(imgstream_fb_write(display, 0, 0, pixel, 1) == 1 &&
	      imgstream_fb_flush(display) == 0,
	      "display framebuffer presents an active image before host teardown");

	delete host;
	CHECK(stream->host() == NULL,
	      "host teardown invalidates registered framebuffer stream host");
	CHECK(stream->present() == -1,
	      "detached framebuffer stream refuses presentation without a host");
	BObolWindowHost replacement;
	CHECK(stream->setHost(&replacement) == 0 &&
	      stream->framebuffer() == stream_fb &&
	      replacement.getFramebufferCount() == 1 &&
	      stream->present() == 0,
	      "stream reattaches after host teardown without reopening storage");
	unsigned char preserved[3] = {0, 0, 0};
	CHECK(stream->read(0, 0, preserved, 1) == 1 &&
	      preserved[0] == pixel[0] && preserved[1] == pixel[1] &&
	      preserved[2] == pixel[2],
	      "host teardown and reattach preserve framebuffer pixels");
	stream->close();
	CHECK(replacement.getFramebufferCount() == 0,
	      "reattached stream close removes replacement-host nodes");
	CHECK(imgstream_fb_flush(display) == 0,
	      "display framebuffer remains usable after host teardown");
	CHECK(imgstream_fb_view(display, 1, 1, 2, 2) == 0,
	      "display framebuffer no longer dispatches to destroyed host");
	imgstream_fb_close(display);
	delete stream;
    }
    return 0;
}

struct DisplaySessionTestProvider {
    BObolWindowHost host;
};

static int display_session_test_polls = 0;

static int
display_session_test_open(bobol_display_endpoint_t *endpoint,
	const imgstream_fb_spec_info_t *spec, size_t UNUSED(width),
	size_t UNUSED(height), const char *UNUSED(title), void **instance,
	void *UNUSED(data))
{
    if (!endpoint || !spec || spec->display != IMGSTREAM_FB_DISPLAY_SWRAST ||
	!instance)
	return 0;
    DisplaySessionTestProvider *provider = new DisplaySessionTestProvider;
    if (!bobol_display_endpoint_host_bind(endpoint, &provider->host, 0)) {
	delete provider;
	return 0;
    }
    *instance = provider;
    return 1;
}

static void
display_session_test_close(void *instance, void *UNUSED(data))
{
    delete static_cast<DisplaySessionTestProvider *>(instance);
}

static int
display_session_test_poll(void *UNUSED(instance), void *UNUSED(data))
{
    display_session_test_polls++;
    return 0;
}

struct DisplaySessionTask {
    imgstream_fb_t *framebuffer = NULL;
    std::thread::id worker_thread;
};

static int
display_session_test_task(void *data)
{
    DisplaySessionTask *task = static_cast<DisplaySessionTask *>(data);
    if (!task || !task->framebuffer)
	return -1;
    task->worker_thread = std::this_thread::get_id();
    const unsigned char pixel[3] = {51, 34, 17};
    if (imgstream_fb_write(task->framebuffer, 1, 1, pixel, 1) != 1 ||
	imgstream_fb_flush(task->framebuffer) != 0)
	return -1;
    std::this_thread::sleep_for(std::chrono::milliseconds(35));
    return 37;
}

static int
test_display_session_contract(void)
{
    static const bobol_display_provider_t provider = {
	BOBOL_DISPLAY_PROVIDER_ABI_VERSION,
	sizeof(bobol_display_provider_t),
	"display-session-test",
	1,
	NULL,
	display_session_test_open,
	display_session_test_close,
	display_session_test_poll,
	NULL
    };
    CHECK(bobol_display_provider_register(&provider),
	  "display session accepts a toolkit-neutral provider");
    CHECK(bobol_display_provider_register(&provider),
	  "display session provider registration is idempotent");

    bobol_display_session_t *session = bobol_display_session_open(
	"/dev/swrast", 4, 3, "display session test");
    CHECK(session != NULL, "display session opens a legacy software target");
    imgstream_fb_t *fb = bobol_display_session_framebuffer(session);
    CHECK(fb != NULL && imgstream_fb_width(fb) == 4 &&
	imgstream_fb_height(fb) == 3,
	"display session returns its image-stream framebuffer");
    const unsigned char pixel[3] = {17, 34, 51};
    CHECK(imgstream_fb_write(fb, 0, 0, pixel, 1) == 1 &&
	imgstream_fb_flush(fb) == 0,
	"display session publishes framebuffer writes through Obol");
    CHECK(imgstream_fb_poll(fb) == 0 && display_session_test_polls > 0,
	"display session polling delegates native event processing to the provider");
    DisplaySessionTask task;
    task.framebuffer = fb;
    display_session_test_polls = 0;
    const std::thread::id owner_thread = std::this_thread::get_id();
    CHECK(bobol_display_session_run(session, display_session_test_task,
	  &task) == 37 && task.worker_thread != owner_thread &&
	  display_session_test_polls > 0,
	"display session keeps its owner thread polling while a worker publishes pixels");
    bobol_display_session_close(session);
    return 0;
}

int
main(int ac, char **av)
{
    bu_setprogname(av[0]);
    (void)ac;
    (void)av;

    bobol_init(NULL);

    if (test_input_context())
	return 1;
    if (test_progressive_status_contract())
	return 1;
    if (test_window_host_contract())
	return 1;
    if (test_display_endpoint_contract())
	return 1;
    if (test_display_endpoint_input())
	return 1;
    if (test_host_factory_contract())
	return 1;
    if (test_host_factory_simultaneous_endpoints())
	return 1;
    if (test_host_factory_reattach())
	return 1;
    if (test_frame_request_callback_teardown())
	return 1;
    if (test_factory_endpoint_close_during_frame_request())
	return 1;
    if (test_imgstream_display_host_bridge())
	return 1;
    if (test_framebuffer_stream_helper())
	return 1;
    if (test_framebuffer_host_teardown())
	return 1;
    if (test_display_session_contract())
	return 1;

    return 0;
}
