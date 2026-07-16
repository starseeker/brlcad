/*                       I N P U T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/input.h"
#include "brlobol/display_endpoint.h"
#include "bv.h"

#include <limits.h>

BRLObolInputEvent::BRLObolInputEvent(void) :
    type(BRLOBOL_INPUT_NONE),
    timestamp(0),
    x(0),
    y(0),
    dx(0),
    dy(0),
    wheelDelta(0),
    button(BRLOBOL_INPUT_ANY),
    buttons(0),
    key(BRLOBOL_INPUT_ANY),
    modifiers(BRLOBOL_INPUT_MOD_NONE),
    width(0),
    height(0)
{
}

BRLObolInputBinding::BRLObolInputBinding(void) :
    eventType(BRLOBOL_INPUT_NONE),
    key(BRLOBOL_INPUT_ANY),
    button(BRLOBOL_INPUT_ANY),
    requiredModifiers(BRLOBOL_INPUT_MOD_NONE),
    forbiddenModifiers(BRLOBOL_INPUT_MOD_NONE),
    priority(0),
    action(BRLOBOL_ACTION_NONE)
{
}

BRLObolInputBinding::BRLObolInputBinding(
    BRLObolInputEventType nextEventType, int nextKey, int nextButton,
    unsigned int nextRequiredModifiers, unsigned int nextForbiddenModifiers,
    int nextPriority, BRLObolInputAction nextAction) :
    eventType(nextEventType),
    key(nextKey),
    button(nextButton),
    requiredModifiers(nextRequiredModifiers),
    forbiddenModifiers(nextForbiddenModifiers),
    priority(nextPriority),
    action(nextAction)
{
}

BRLObolInputContext::BRLObolInputContext(void) :
    handler(NULL),
    handlerData(NULL),
    layerHandler(NULL),
    layerHandlerData(NULL),
    layerOwner(NULL)
{
}

BRLObolInputContext::~BRLObolInputContext(void)
{
}

void
BRLObolInputContext::clear(void)
{
    this->bindings.clear();
    this->layerBindings.clear();
    this->handler = NULL;
    this->handlerData = NULL;
    this->layerHandler = NULL;
    this->layerHandlerData = NULL;
    this->layerOwner = NULL;
}

void
BRLObolInputContext::setProfile(const BRLObolInputProfile *profile)
{
    this->setBindings(profile ? profile->bindings : NULL,
	profile ? profile->bindingCount : 0);
}

void
BRLObolInputContext::setBindings(const BRLObolInputBinding *nextBindings,
	size_t count)
{
    this->bindings.clear();
    if (nextBindings && count)
	this->bindings.assign(nextBindings, nextBindings + count);
}

void
BRLObolInputContext::setActionHandler(BRLObolInputActionHandler nextHandler,
	void *userData)
{
    this->handler = nextHandler;
    this->handlerData = userData;
}

int
BRLObolInputContext::clearActionHandlerIf(
	BRLObolInputActionHandler expectedHandler, void *expectedUserData)
{
    if (this->handler != expectedHandler ||
	this->handlerData != expectedUserData)
	return 0;
    this->handler = NULL;
    this->handlerData = NULL;
    return 1;
}

int
BRLObolInputContext::setActionLayer(const BRLObolInputActionLayer *layer,
	void *owner, void *userData)
{
    if (!layer || !layer->handler || !layer->bindings ||
	!layer->bindingCount || !owner ||
	(this->layerOwner && this->layerOwner != owner))
	return 0;

    this->layerBindings.assign(layer->bindings,
	layer->bindings + layer->bindingCount);
    this->layerHandler = layer->handler;
    this->layerHandlerData = userData;
    this->layerOwner = owner;
    return 1;
}

int
BRLObolInputContext::clearActionLayerIf(void *owner)
{
    if (!owner || this->layerOwner != owner)
	return 0;
    this->layerBindings.clear();
    this->layerHandler = NULL;
    this->layerHandlerData = NULL;
    this->layerOwner = NULL;
    return 1;
}

static unsigned int
modifier_count(unsigned int modifiers)
{
    unsigned int count = 0;
    while (modifiers) {
	count += modifiers & 1u;
	modifiers >>= 1;
    }
    return count;
}

static int
binding_matches(const BRLObolInputBinding &binding,
	const BRLObolInputEvent &event)
{
    if (binding.eventType != event.type)
	return 0;
    if (binding.key != BRLOBOL_INPUT_ANY && binding.key != event.key)
	return 0;
    if (binding.button != BRLOBOL_INPUT_ANY &&
	binding.button != event.button)
	return 0;
    if ((event.modifiers & binding.requiredModifiers) !=
	binding.requiredModifiers)
	return 0;
    return (event.modifiers & binding.forbiddenModifiers) == 0 ? 1 : 0;
}

int
BRLObolInputContext::dispatch(const BRLObolInputEvent *event) const
{
    if (!event)
	return -1;

    const BRLObolInputBinding *baseWinner = NULL;
    const BRLObolInputBinding *layerWinner = NULL;
    int baseScore = INT_MIN;
    int layerScore = INT_MIN;

    const auto consider = [event](const BRLObolInputBinding &binding,
	const BRLObolInputBinding **winner, int *winnerScore) {
	if (!binding_matches(binding, *event))
	    return;
	const int score = binding.priority * 32 +
	    static_cast<int>(modifier_count(binding.requiredModifiers));
	if (!*winner || score > *winnerScore) {
	    *winner = &binding;
	    *winnerScore = score;
	}
    };
    for (const BRLObolInputBinding &binding : this->bindings)
	consider(binding, &baseWinner, &baseScore);
    for (const BRLObolInputBinding &binding : this->layerBindings)
	consider(binding, &layerWinner, &layerScore);
    if (!baseWinner && !layerWinner)
	return 0;

    /* A handler may replace its owner-scoped layer while beginning or ending
     * a gesture.  Copy both candidate actions and handlers before invoking it
     * so replacement cannot invalidate a binding needed for fallthrough. */
    const BRLObolInputAction baseAction = baseWinner ? baseWinner->action :
	BRLOBOL_ACTION_NONE;
    const BRLObolInputAction layerAction = layerWinner ?
	layerWinner->action : BRLOBOL_ACTION_NONE;
    BRLObolInputActionHandler baseHandler = this->handler;
    void *baseHandlerData = this->handlerData;
    BRLObolInputActionHandler appHandler = this->layerHandler;
    void *appHandlerData = this->layerHandlerData;

    if (layerWinner && layerScore >= baseScore &&
	layerAction != BRLOBOL_ACTION_NONE) {
	const int result = appHandler ? appHandler(appHandlerData, layerAction,
	    event) : BRLOBOL_INPUT_RESULT_UNHANDLED;
	if (result != BRLOBOL_INPUT_RESULT_UNHANDLED)
	    return result;
    }
    return baseWinner && baseAction != BRLOBOL_ACTION_NONE && baseHandler ?
	baseHandler(baseHandlerData, baseAction, event) :
	BRLOBOL_INPUT_RESULT_UNHANDLED;
}

int
BRLObolInputContext::hasAction(BRLObolInputAction action) const
{
    for (const BRLObolInputBinding &binding : this->bindings) {
	if (binding.action == action)
	    return 1;
    }
    for (const BRLObolInputBinding &binding : this->layerBindings) {
	if (binding.action == action)
	    return 1;
    }
    return 0;
}

size_t
BRLObolInputContext::bindingCount(void) const
{
    return this->bindings.size() + this->layerBindings.size();
}

extern "C" const BRLObolInputProfile *
brlobol_input_default_view_profile(void)
{
    /* View shortcuts are intentionally unmodified.  Shift has separate MGED
     * bindings such as accept/rear, and control/alt/meta belong to the host.
     */
    static const unsigned int plain_modifiers = BRLOBOL_INPUT_MOD_SHIFT |
	BRLOBOL_INPUT_MOD_CONTROL | BRLOBOL_INPUT_MOD_ALT |
	BRLOBOL_INPUT_MOD_META;
    static const BRLObolInputBinding bindings[] = {
	{BRLOBOL_INPUT_KEY_PRESS, 'A', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_TOGGLE_ADC},
	{BRLOBOL_INPUT_KEY_PRESS, 'M', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_TOGGLE_MODEL_AXES},
	{BRLOBOL_INPUT_KEY_PRESS, 'V', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_TOGGLE_VIEW_AXES},
	{BRLOBOL_INPUT_KEY_PRESS, '2', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_2},
	{BRLOBOL_INPUT_KEY_PRESS, '3', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_3},
	{BRLOBOL_INPUT_KEY_PRESS, '4', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_4},
	{BRLOBOL_INPUT_KEY_PRESS, '5', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_5},
	{BRLOBOL_INPUT_KEY_PRESS, '6', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_6},
	{BRLOBOL_INPUT_KEY_PRESS, '7', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_7},
	{BRLOBOL_INPUT_KEY_PRESS, 'F', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_FRONT},
	{BRLOBOL_INPUT_KEY_PRESS, 'T', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_TOP},
	{BRLOBOL_INPUT_KEY_PRESS, 'B', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_BOTTOM},
	{BRLOBOL_INPUT_KEY_PRESS, 'L', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_LEFT},
	{BRLOBOL_INPUT_KEY_PRESS, 'R', BRLOBOL_INPUT_ANY,
	 BRLOBOL_INPUT_MOD_SHIFT,
	 BRLOBOL_INPUT_MOD_CONTROL | BRLOBOL_INPUT_MOD_ALT |
	 BRLOBOL_INPUT_MOD_META, 1, BRLOBOL_ACTION_VIEW_REAR},
	{BRLOBOL_INPUT_KEY_PRESS, 'R', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_RIGHT},
	{BRLOBOL_INPUT_POINTER_PRESS, BRLOBOL_INPUT_ANY, 0,
	 BRLOBOL_INPUT_MOD_SHIFT,
	 BRLOBOL_INPUT_MOD_CONTROL | BRLOBOL_INPUT_MOD_ALT |
	 BRLOBOL_INPUT_MOD_META, 1, BRLOBOL_ACTION_VIEW_PAN_BEGIN},
	{BRLOBOL_INPUT_POINTER_RELEASE, BRLOBOL_INPUT_ANY, 0,
	 BRLOBOL_INPUT_MOD_SHIFT,
	 BRLOBOL_INPUT_MOD_CONTROL | BRLOBOL_INPUT_MOD_ALT |
	 BRLOBOL_INPUT_MOD_META, 1, BRLOBOL_ACTION_VIEW_PAN_END},
	{BRLOBOL_INPUT_POINTER_RELEASE, BRLOBOL_INPUT_ANY, 0, 0,
	 BRLOBOL_INPUT_MOD_SHIFT | BRLOBOL_INPUT_MOD_CONTROL |
	 BRLOBOL_INPUT_MOD_ALT | BRLOBOL_INPUT_MOD_META, 0,
	 BRLOBOL_ACTION_VIEW_PRIMARY_RELEASE},
	{BRLOBOL_INPUT_POINTER_RELEASE, BRLOBOL_INPUT_ANY, 1, 0,
	 BRLOBOL_INPUT_MOD_SHIFT | BRLOBOL_INPUT_MOD_CONTROL |
	 BRLOBOL_INPUT_MOD_ALT | BRLOBOL_INPUT_MOD_META, 0,
	 BRLOBOL_ACTION_VIEW_CENTER_RELEASE},
	{BRLOBOL_INPUT_POINTER_RELEASE, BRLOBOL_INPUT_ANY, 2, 0,
	 BRLOBOL_INPUT_MOD_SHIFT | BRLOBOL_INPUT_MOD_CONTROL |
	 BRLOBOL_INPUT_MOD_ALT | BRLOBOL_INPUT_MOD_META, 0,
	 BRLOBOL_ACTION_VIEW_SECONDARY_RELEASE},
	{BRLOBOL_INPUT_POINTER_MOTION, BRLOBOL_INPUT_ANY, 0,
	 BRLOBOL_INPUT_MOD_CONTROL | BRLOBOL_INPUT_MOD_SHIFT, 0, 2,
	 BRLOBOL_ACTION_VIEW_ZOOM},
	{BRLOBOL_INPUT_POINTER_MOTION, BRLOBOL_INPUT_ANY, 0,
	 BRLOBOL_INPUT_MOD_CONTROL, BRLOBOL_INPUT_MOD_SHIFT, 1,
	 BRLOBOL_ACTION_VIEW_ROTATE},
	{BRLOBOL_INPUT_POINTER_MOTION, BRLOBOL_INPUT_ANY, 0,
	 BRLOBOL_INPUT_MOD_SHIFT, BRLOBOL_INPUT_MOD_CONTROL, 1,
	 BRLOBOL_ACTION_VIEW_PAN},
	{BRLOBOL_INPUT_POINTER_MOTION, BRLOBOL_INPUT_ANY, 0, 0,
	 BRLOBOL_INPUT_MOD_CONTROL | BRLOBOL_INPUT_MOD_SHIFT, 0,
	 BRLOBOL_ACTION_VIEW_ADJUST},
	{BRLOBOL_INPUT_WHEEL, BRLOBOL_INPUT_ANY, BRLOBOL_INPUT_ANY, 0, 0, 0,
	 BRLOBOL_ACTION_VIEW_ZOOM}
    };
    static const BRLObolInputProfile profile = {
	"brlcad-default-view", bindings,
	sizeof(bindings) / sizeof(bindings[0])
    };
    return &profile;
}

extern "C" const BRLObolInputProfile *
brlobol_input_keyboard_view_profile(void)
{
    /* Keep pointer gestures application-owned when their mode depends on
     * editing or selection state outside the common view controller. */
    static const unsigned int plain_modifiers = BRLOBOL_INPUT_MOD_SHIFT |
	BRLOBOL_INPUT_MOD_CONTROL | BRLOBOL_INPUT_MOD_ALT |
	BRLOBOL_INPUT_MOD_META;
    static const BRLObolInputBinding bindings[] = {
	{BRLOBOL_INPUT_KEY_PRESS, 'A', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_TOGGLE_ADC},
	{BRLOBOL_INPUT_KEY_PRESS, 'M', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_TOGGLE_MODEL_AXES},
	{BRLOBOL_INPUT_KEY_PRESS, 'V', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_TOGGLE_VIEW_AXES},
	{BRLOBOL_INPUT_KEY_PRESS, '2', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_2},
	{BRLOBOL_INPUT_KEY_PRESS, '3', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_3},
	{BRLOBOL_INPUT_KEY_PRESS, '4', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_4},
	{BRLOBOL_INPUT_KEY_PRESS, '5', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_5},
	{BRLOBOL_INPUT_KEY_PRESS, '6', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_6},
	{BRLOBOL_INPUT_KEY_PRESS, '7', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_7},
	{BRLOBOL_INPUT_KEY_PRESS, 'F', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_FRONT},
	{BRLOBOL_INPUT_KEY_PRESS, 'T', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_TOP},
	{BRLOBOL_INPUT_KEY_PRESS, 'B', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_BOTTOM},
	{BRLOBOL_INPUT_KEY_PRESS, 'L', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_LEFT},
	{BRLOBOL_INPUT_KEY_PRESS, 'R', BRLOBOL_INPUT_ANY,
	 BRLOBOL_INPUT_MOD_SHIFT,
	 BRLOBOL_INPUT_MOD_CONTROL | BRLOBOL_INPUT_MOD_ALT |
	 BRLOBOL_INPUT_MOD_META, 1, BRLOBOL_ACTION_VIEW_REAR},
	{BRLOBOL_INPUT_KEY_PRESS, 'R', BRLOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BRLOBOL_ACTION_VIEW_RIGHT}
    };
    static const BRLObolInputProfile profile = {
	"brlcad-keyboard-view", bindings,
	sizeof(bindings) / sizeof(bindings[0])
    };
    return &profile;
}

extern "C" int
brlobol_input_view_orientation_apply(void *view_ctx,
	BRLObolInputAction action)
{
    struct bv_context *context = static_cast<struct bv_context *>(view_ctx);
    struct bv *view = bv_context_view(context);
    if (!view)
	return 0;

    fastf_t azimuth = 0.0;
    fastf_t elevation = 0.0;
    switch (action) {
	case BRLOBOL_ACTION_VIEW_2:
	    azimuth = 35.0;
	    elevation = -25.0;
	    break;
	case BRLOBOL_ACTION_VIEW_3:
	    azimuth = 35.0;
	    elevation = 25.0;
	    break;
	case BRLOBOL_ACTION_VIEW_4:
	    azimuth = 45.0;
	    elevation = 45.0;
	    break;
	case BRLOBOL_ACTION_VIEW_5:
	    azimuth = 145.0;
	    elevation = 25.0;
	    break;
	case BRLOBOL_ACTION_VIEW_6:
	    azimuth = 215.0;
	    elevation = 25.0;
	    break;
	case BRLOBOL_ACTION_VIEW_7:
	    azimuth = 325.0;
	    elevation = 25.0;
	    break;
	case BRLOBOL_ACTION_VIEW_TOP:
	    azimuth = 270.0;
	    elevation = 90.0;
	    break;
	case BRLOBOL_ACTION_VIEW_BOTTOM:
	    azimuth = 270.0;
	    elevation = -90.0;
	    break;
	case BRLOBOL_ACTION_VIEW_LEFT:
	    azimuth = 90.0;
	    break;
	case BRLOBOL_ACTION_VIEW_REAR:
	    azimuth = 180.0;
	    break;
	case BRLOBOL_ACTION_VIEW_RIGHT:
	    azimuth = 270.0;
	    break;
	case BRLOBOL_ACTION_VIEW_FRONT:
	    break;
	default:
	    return 0;
    }

    vect_t aet;
    VSET(aet, azimuth, elevation, 0.0);
    if (!bv_aet_set(view, aet))
	return 0;
    return bv_context_update(context, BV_CONTEXT_CHANGED_VIEW) ? 1 : 0;
}

static const char *
brlobol_input_faceplate_visibility_property(BRLObolInputAction action)
{
    switch (action) {
	case BRLOBOL_ACTION_TOGGLE_ADC:
	    return "view.faceplate.adc.visible";
	case BRLOBOL_ACTION_TOGGLE_MODEL_AXES:
	    return "view.faceplate.model_axes.visible";
	case BRLOBOL_ACTION_TOGGLE_VIEW_AXES:
	    return "view.faceplate.view_axes.visible";
	default:
	    return NULL;
    }
}

extern "C" int
brlobol_display_endpoint_input_faceplate_toggle_apply(
	brlobol_display_endpoint_t *endpoint, void *view_ctx,
	BRLObolInputAction action, int *visible)
{
    const char *property_name =
	brlobol_input_faceplate_visibility_property(action);
    if (!property_name)
	return 0;

    if (endpoint) {
	struct brlobol_endpoint_property_value value =
	    BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	const int get_result = brlobol_display_endpoint_property_get(endpoint,
	    property_name, &value);
	if (get_result == BRLOBOL_ENDPOINT_PROPERTY_OK) {
	    if (value.type != BRLOBOL_ENDPOINT_PROPERTY_BOOL)
		return 0;
	    value.bool_value = value.bool_value ? 0 : 1;
	    const int set_result = brlobol_display_endpoint_property_set(endpoint,
		property_name, &value);
	    if (set_result != BRLOBOL_ENDPOINT_PROPERTY_OK)
		return 0;
	    if (visible)
		*visible = value.bool_value;
	    return 1;
	}
	/* A standalone controller has no GED property owner.  It retains the
	 * same state locally rather than inventing an application renderer path. */
	if (get_result != BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED)
	    return 0;
    }

    if (!view_ctx)
	return 0;
    struct bv_context *context = static_cast<struct bv_context *>(view_ctx);
    struct bv *view = bv_context_view(context);
    if (!view)
	return 0;

    if (action == BRLOBOL_ACTION_TOGGLE_ADC) {
	struct bv_adc_state state = {};
	if (!bv_adc_state_get(&state, view))
	    return 0;
	state.draw = state.draw ? 0 : 1;
	if (!bv_adc_state_set(view, &state))
	    return 0;
	if (visible)
	    *visible = state.draw ? 1 : 0;
	return 1;
    }

    struct bv_axes_state state = {};
    const int model_axes = action == BRLOBOL_ACTION_TOGGLE_MODEL_AXES;
    if (!(model_axes ? bv_model_axes_state_get(&state, view) :
	bv_view_axes_state_get(&state, view)))
	return 0;
    state.draw = state.draw ? 0 : 1;
    if (!(model_axes ? bv_model_axes_state_set(view, &state) :
	bv_view_axes_state_set(view, &state)))
	return 0;
    if (visible)
	*visible = state.draw ? 1 : 0;
    return 1;
}

const BRLObolInputProfile &
BRLObolInputContext::defaultViewProfile(void)
{
    return *brlobol_input_default_view_profile();
}
