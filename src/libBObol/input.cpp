/*                       I N P U T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BInput.h"
#include "BObol/BDisplayEndpoint.h"
#include "identity_counter_private.h"
#include "bv.h"

#include <algorithm>
#include <limits.h>

BObolInputEvent::BObolInputEvent(void) :
    type(BOBOL_INPUT_NONE),
    timestamp(0),
    x(0),
    y(0),
    dx(0),
    dy(0),
    wheelDelta(0),
    button(BOBOL_INPUT_ANY),
    buttons(0),
    key(BOBOL_INPUT_ANY),
    modifiers(BOBOL_INPUT_MOD_NONE),
    width(0),
    height(0)
{
}

BObolInputBinding::BObolInputBinding(void) :
    eventType(BOBOL_INPUT_NONE),
    key(BOBOL_INPUT_ANY),
    button(BOBOL_INPUT_ANY),
    requiredModifiers(BOBOL_INPUT_MOD_NONE),
    forbiddenModifiers(BOBOL_INPUT_MOD_NONE),
    priority(0),
    action(BOBOL_ACTION_NONE)
{
}

BObolInputBinding::BObolInputBinding(
    BObolInputEventType nextEventType, int nextKey, int nextButton,
    unsigned int nextRequiredModifiers, unsigned int nextForbiddenModifiers,
    int nextPriority, BObolInputAction nextAction) :
    eventType(nextEventType),
    key(nextKey),
    button(nextButton),
    requiredModifiers(nextRequiredModifiers),
    forbiddenModifiers(nextForbiddenModifiers),
    priority(nextPriority),
    action(nextAction)
{
}

BObolInputContext::BObolInputContext(void) :
    handler(NULL),
    handlerData(NULL),
    nextLayerOrder(1)
{
}

BObolInputContext::~BObolInputContext(void)
{
}

void
BObolInputContext::clear(void)
{
    this->bindings.clear();
    this->actionLayers.clear();
    this->handler = NULL;
    this->handlerData = NULL;
    this->nextLayerOrder = 1;
}

void
BObolInputContext::setProfile(const BObolInputProfile *profile)
{
    this->setBindings(profile ? profile->bindings : NULL,
	profile ? profile->bindingCount : 0);
}

void
BObolInputContext::setBindings(const BObolInputBinding *nextBindings,
	size_t count)
{
    this->bindings.clear();
    if (nextBindings && count)
	this->bindings.assign(nextBindings, nextBindings + count);
}

void
BObolInputContext::setActionHandler(BObolInputActionHandler nextHandler,
	void *userData)
{
    this->handler = nextHandler;
    this->handlerData = userData;
}

int
BObolInputContext::clearActionHandlerIf(
	BObolInputActionHandler expectedHandler, void *expectedUserData)
{
    if (this->handler != expectedHandler ||
	this->handlerData != expectedUserData)
	return 0;
    this->handler = NULL;
    this->handlerData = NULL;
    return 1;
}

int
BObolInputContext::setActionLayer(const BObolInputActionLayer *layer,
	void *owner, void *userData)
{
    if (!layer || !layer->handler || !layer->bindings ||
	!layer->bindingCount || !owner)
	return 0;

    ActionLayer *target = NULL;
    for (ActionLayer &entry : this->actionLayers) {
	if (entry.owner == owner) {
	    target = &entry;
	    break;
	}
    }
    if (!target) {
	this->actionLayers.push_back(ActionLayer());
	target = &this->actionLayers.back();
	target->owner = owner;
    }
    target->bindings.assign(layer->bindings,
	layer->bindings + layer->bindingCount);
    target->handler = layer->handler;
    target->handlerData = userData;
    target->order = bobol_nonzero_identity_take(this->nextLayerOrder);
    return 1;
}

int
BObolInputContext::clearActionLayerIf(void *owner)
{
    if (!owner)
	return 0;
    for (std::vector<ActionLayer>::iterator it = this->actionLayers.begin();
	 it != this->actionLayers.end(); ++it) {
	if (it->owner != owner)
	    continue;
	this->actionLayers.erase(it);
	return 1;
    }
    return 0;
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
binding_matches(const BObolInputBinding &binding,
	const BObolInputEvent &event)
{
    if (binding.eventType != event.type)
	return 0;
    if (binding.key != BOBOL_INPUT_ANY && binding.key != event.key)
	return 0;
    if (binding.button != BOBOL_INPUT_ANY &&
	binding.button != event.button)
	return 0;
    if ((event.modifiers & binding.requiredModifiers) !=
	binding.requiredModifiers)
	return 0;
    return (event.modifiers & binding.forbiddenModifiers) == 0 ? 1 : 0;
}

int
BObolInputContext::dispatch(const BObolInputEvent *event) const
{
    if (!event)
	return -1;

    struct Candidate {
	BObolInputAction action;
	BObolInputActionHandler handler;
	void *handlerData;
	int score;
	uint64_t order;
	int layer;
    };
    std::vector<Candidate> candidates;

    const auto consider = [event](const BObolInputBinding &binding,
	const BObolInputBinding **winner, int *winnerScore) {
	if (!binding_matches(binding, *event))
	    return;
	const int score = binding.priority * 32 +
	    static_cast<int>(modifier_count(binding.requiredModifiers));
	if (!*winner || score > *winnerScore) {
	    *winner = &binding;
	    *winnerScore = score;
	}
    };
    const BObolInputBinding *baseWinner = NULL;
    int baseScore = INT_MIN;
    for (const BObolInputBinding &binding : this->bindings)
	consider(binding, &baseWinner, &baseScore);
    if (baseWinner && baseWinner->action != BOBOL_ACTION_NONE &&
	this->handler) {
	Candidate candidate = {baseWinner->action, this->handler,
	    this->handlerData, baseScore, 0, 0};
	candidates.push_back(candidate);
    }

    for (const ActionLayer &layer : this->actionLayers) {
	const BObolInputBinding *winner = NULL;
	int score = INT_MIN;
	for (const BObolInputBinding &binding : layer.bindings)
	    consider(binding, &winner, &score);
	if (!winner || winner->action == BOBOL_ACTION_NONE || !layer.handler)
	    continue;
	Candidate candidate = {winner->action, layer.handler,
	    layer.handlerData, score, layer.order, 1};
	candidates.push_back(candidate);
    }

    std::stable_sort(candidates.begin(), candidates.end(),
	[](const Candidate &a, const Candidate &b) {
	    if (a.score != b.score)
		return a.score > b.score;
	    if (a.layer != b.layer)
		return a.layer > b.layer;
	    return a.order > b.order;
	});

    /* Handlers may replace or remove layers during dispatch.  Candidates are
     * value snapshots, so fallthrough never dereferences mutated storage. */
    for (const Candidate &candidate : candidates) {
	const int result = candidate.handler(candidate.handlerData,
	    candidate.action, event);
	if (result != BOBOL_INPUT_RESULT_UNHANDLED)
	    return result;
    }
    return BOBOL_INPUT_RESULT_UNHANDLED;
}

int
BObolInputContext::hasAction(BObolInputAction action) const
{
    for (const BObolInputBinding &binding : this->bindings) {
	if (binding.action == action)
	    return 1;
    }
    for (const ActionLayer &layer : this->actionLayers) {
	for (const BObolInputBinding &binding : layer.bindings) {
	    if (binding.action == action)
		return 1;
	}
    }
    return 0;
}

size_t
BObolInputContext::bindingCount(void) const
{
    size_t count = this->bindings.size();
    for (const ActionLayer &layer : this->actionLayers)
	count += layer.bindings.size();
    return count;
}

extern "C" const BObolInputProfile *
bobol_input_default_view_profile(void)
{
    /* View shortcuts are intentionally unmodified.  Shift has separate MGED
     * bindings such as accept/rear, and control/alt/meta belong to the host.
     */
    static const unsigned int plain_modifiers = BOBOL_INPUT_MOD_SHIFT |
	BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_ALT |
	BOBOL_INPUT_MOD_META;
    static const BObolInputBinding bindings[] = {
	{BOBOL_INPUT_KEY_PRESS, 'A', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_TOGGLE_ADC},
	{BOBOL_INPUT_KEY_PRESS, 'M', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_TOGGLE_MODEL_AXES},
	{BOBOL_INPUT_KEY_PRESS, 'V', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_TOGGLE_VIEW_AXES},
	{BOBOL_INPUT_KEY_PRESS, '2', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_2},
	{BOBOL_INPUT_KEY_PRESS, '3', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_3},
	{BOBOL_INPUT_KEY_PRESS, '4', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_4},
	{BOBOL_INPUT_KEY_PRESS, '5', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_5},
	{BOBOL_INPUT_KEY_PRESS, '6', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_6},
	{BOBOL_INPUT_KEY_PRESS, '7', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_7},
	{BOBOL_INPUT_KEY_PRESS, 'F', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_FRONT},
	{BOBOL_INPUT_KEY_PRESS, 'T', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_TOP},
	{BOBOL_INPUT_KEY_PRESS, 'B', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_BOTTOM},
	{BOBOL_INPUT_KEY_PRESS, 'L', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_LEFT},
	{BOBOL_INPUT_KEY_PRESS, 'R', BOBOL_INPUT_ANY,
	 BOBOL_INPUT_MOD_SHIFT,
	 BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_ALT |
	 BOBOL_INPUT_MOD_META, 1, BOBOL_ACTION_VIEW_REAR},
	{BOBOL_INPUT_KEY_PRESS, 'R', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_RIGHT},
	{BOBOL_INPUT_POINTER_PRESS, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_SHIFT,
	 BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_ALT |
	 BOBOL_INPUT_MOD_META, 1, BOBOL_ACTION_VIEW_PAN_BEGIN},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_SHIFT,
	 BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_ALT |
	 BOBOL_INPUT_MOD_META, 1, BOBOL_ACTION_VIEW_PAN_END},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, 0, 0,
	 BOBOL_INPUT_MOD_SHIFT | BOBOL_INPUT_MOD_CONTROL |
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META, 0,
	 BOBOL_ACTION_VIEW_PRIMARY_RELEASE},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, 1, 0,
	 BOBOL_INPUT_MOD_SHIFT | BOBOL_INPUT_MOD_CONTROL |
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META, 0,
	 BOBOL_ACTION_VIEW_CENTER_RELEASE},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, 2, 0,
	 BOBOL_INPUT_MOD_SHIFT | BOBOL_INPUT_MOD_CONTROL |
	 BOBOL_INPUT_MOD_ALT | BOBOL_INPUT_MOD_META, 0,
	 BOBOL_ACTION_VIEW_SECONDARY_RELEASE},
	{BOBOL_INPUT_POINTER_MOTION, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_SHIFT, 0, 2,
	 BOBOL_ACTION_VIEW_ZOOM},
	{BOBOL_INPUT_POINTER_MOTION, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_CONTROL, BOBOL_INPUT_MOD_SHIFT, 1,
	 BOBOL_ACTION_VIEW_ROTATE},
	{BOBOL_INPUT_POINTER_MOTION, BOBOL_INPUT_ANY, 0,
	 BOBOL_INPUT_MOD_SHIFT, BOBOL_INPUT_MOD_CONTROL, 1,
	 BOBOL_ACTION_VIEW_PAN},
	{BOBOL_INPUT_POINTER_MOTION, BOBOL_INPUT_ANY, 0, 0,
	 BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_SHIFT, 0,
	 BOBOL_ACTION_VIEW_ADJUST},
	{BOBOL_INPUT_WHEEL, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY, 0, 0, 0,
	 BOBOL_ACTION_VIEW_ZOOM}
    };
    static const BObolInputProfile profile = {
	"brlcad-default-view", bindings,
	sizeof(bindings) / sizeof(bindings[0])
    };
    return &profile;
}

extern "C" const BObolInputProfile *
bobol_input_keyboard_view_profile(void)
{
    /* Keep pointer gestures application-owned when their mode depends on
     * editing or selection state outside the common view controller. */
    static const unsigned int plain_modifiers = BOBOL_INPUT_MOD_SHIFT |
	BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_ALT |
	BOBOL_INPUT_MOD_META;
    static const BObolInputBinding bindings[] = {
	{BOBOL_INPUT_KEY_PRESS, 'A', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_TOGGLE_ADC},
	{BOBOL_INPUT_KEY_PRESS, 'M', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_TOGGLE_MODEL_AXES},
	{BOBOL_INPUT_KEY_PRESS, 'V', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_TOGGLE_VIEW_AXES},
	{BOBOL_INPUT_KEY_PRESS, '2', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_2},
	{BOBOL_INPUT_KEY_PRESS, '3', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_3},
	{BOBOL_INPUT_KEY_PRESS, '4', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_4},
	{BOBOL_INPUT_KEY_PRESS, '5', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_5},
	{BOBOL_INPUT_KEY_PRESS, '6', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_6},
	{BOBOL_INPUT_KEY_PRESS, '7', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_7},
	{BOBOL_INPUT_KEY_PRESS, 'F', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_FRONT},
	{BOBOL_INPUT_KEY_PRESS, 'T', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_TOP},
	{BOBOL_INPUT_KEY_PRESS, 'B', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_BOTTOM},
	{BOBOL_INPUT_KEY_PRESS, 'L', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_LEFT},
	{BOBOL_INPUT_KEY_PRESS, 'R', BOBOL_INPUT_ANY,
	 BOBOL_INPUT_MOD_SHIFT,
	 BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_ALT |
	 BOBOL_INPUT_MOD_META, 1, BOBOL_ACTION_VIEW_REAR},
	{BOBOL_INPUT_KEY_PRESS, 'R', BOBOL_INPUT_ANY, 0, plain_modifiers, 0,
	 BOBOL_ACTION_VIEW_RIGHT}
    };
    static const BObolInputProfile profile = {
	"brlcad-keyboard-view", bindings,
	sizeof(bindings) / sizeof(bindings[0])
    };
    return &profile;
}

extern "C" int
bobol_input_view_orientation_apply(void *view_ctx,
	BObolInputAction action)
{
    struct bv_context *context = static_cast<struct bv_context *>(view_ctx);
    struct bv *view = bv_context_view(context);
    if (!view)
	return 0;

    fastf_t azimuth = 0.0;
    fastf_t elevation = 0.0;
    switch (action) {
	case BOBOL_ACTION_VIEW_2:
	    azimuth = 35.0;
	    elevation = -25.0;
	    break;
	case BOBOL_ACTION_VIEW_3:
	    azimuth = 35.0;
	    elevation = 25.0;
	    break;
	case BOBOL_ACTION_VIEW_4:
	    azimuth = 45.0;
	    elevation = 45.0;
	    break;
	case BOBOL_ACTION_VIEW_5:
	    azimuth = 145.0;
	    elevation = 25.0;
	    break;
	case BOBOL_ACTION_VIEW_6:
	    azimuth = 215.0;
	    elevation = 25.0;
	    break;
	case BOBOL_ACTION_VIEW_7:
	    azimuth = 325.0;
	    elevation = 25.0;
	    break;
	case BOBOL_ACTION_VIEW_TOP:
	    azimuth = 270.0;
	    elevation = 90.0;
	    break;
	case BOBOL_ACTION_VIEW_BOTTOM:
	    azimuth = 270.0;
	    elevation = -90.0;
	    break;
	case BOBOL_ACTION_VIEW_LEFT:
	    azimuth = 90.0;
	    break;
	case BOBOL_ACTION_VIEW_REAR:
	    azimuth = 180.0;
	    break;
	case BOBOL_ACTION_VIEW_RIGHT:
	    azimuth = 270.0;
	    break;
	case BOBOL_ACTION_VIEW_FRONT:
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
bobol_input_faceplate_visibility_property(BObolInputAction action)
{
    switch (action) {
	case BOBOL_ACTION_TOGGLE_ADC:
	    return "view.faceplate.adc.visible";
	case BOBOL_ACTION_TOGGLE_MODEL_AXES:
	    return "view.faceplate.model_axes.visible";
	case BOBOL_ACTION_TOGGLE_VIEW_AXES:
	    return "view.faceplate.view_axes.visible";
	default:
	    return NULL;
    }
}

extern "C" int
bobol_display_endpoint_input_faceplate_toggle_apply(
	bobol_display_endpoint_t *endpoint, void *view_ctx,
	BObolInputAction action, int *visible)
{
    const char *property_name =
	bobol_input_faceplate_visibility_property(action);
    if (!property_name)
	return 0;

    if (endpoint) {
	struct bv_display_property_value value =
	    BV_DISPLAY_PROPERTY_VALUE_INIT;
	const int get_result = bobol_display_endpoint_property_get(endpoint,
	    property_name, &value);
	if (get_result == BV_DISPLAY_PROPERTY_OK) {
	    if (value.type != BV_DISPLAY_PROPERTY_BOOL)
		return 0;
	    value.bool_value = value.bool_value ? 0 : 1;
	    const int set_result = bobol_display_endpoint_property_set(endpoint,
		property_name, &value);
	    if (set_result != BV_DISPLAY_PROPERTY_OK)
		return 0;
	    if (visible)
		*visible = value.bool_value;
	    return 1;
	}
	/* A standalone controller has no GED property owner.  It retains the
	 * same state locally rather than inventing an application renderer path. */
	if (get_result != BV_DISPLAY_PROPERTY_UNSUPPORTED)
	    return 0;
    }

    if (!view_ctx)
	return 0;
    struct bv_context *context = static_cast<struct bv_context *>(view_ctx);
    struct bv *view = bv_context_view(context);
    if (!view)
	return 0;

    if (action == BOBOL_ACTION_TOGGLE_ADC) {
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
    const int model_axes = action == BOBOL_ACTION_TOGGLE_MODEL_AXES;
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

const BObolInputProfile &
BObolInputContext::defaultViewProfile(void)
{
    return *bobol_input_default_view_profile();
}
