/*                       I N P U T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/input.h"

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
    handlerData(NULL)
{
}

BRLObolInputContext::~BRLObolInputContext(void)
{
}

void
BRLObolInputContext::clear(void)
{
    this->bindings.clear();
    this->handler = NULL;
    this->handlerData = NULL;
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

    const BRLObolInputBinding *winner = NULL;
    int winnerScore = INT_MIN;
    for (const BRLObolInputBinding &binding : this->bindings) {
	if (!binding_matches(binding, *event))
	    continue;
	const int score = binding.priority * 32 +
	    static_cast<int>(modifier_count(binding.requiredModifiers));
	if (!winner || score > winnerScore) {
	    winner = &binding;
	    winnerScore = score;
	}
    }
    if (!winner || winner->action == BRLOBOL_ACTION_NONE)
	return 0;
    return this->handler ? this->handler(this->handlerData, winner->action,
	event) : 0;
}

int
BRLObolInputContext::hasAction(BRLObolInputAction action) const
{
    for (const BRLObolInputBinding &binding : this->bindings) {
	if (binding.action == action)
	    return 1;
    }
    return 0;
}

size_t
BRLObolInputContext::bindingCount(void) const
{
    return this->bindings.size();
}

extern "C" const BRLObolInputProfile *
brlobol_input_default_view_profile(void)
{
    static const BRLObolInputBinding bindings[] = {
	{BRLOBOL_INPUT_KEY_PRESS, 'A', BRLOBOL_INPUT_ANY, 0, 0, 0,
	 BRLOBOL_ACTION_TOGGLE_ADC},
	{BRLOBOL_INPUT_KEY_PRESS, 'M', BRLOBOL_INPUT_ANY, 0, 0, 0,
	 BRLOBOL_ACTION_TOGGLE_MODEL_AXES},
	{BRLOBOL_INPUT_KEY_PRESS, 'V', BRLOBOL_INPUT_ANY, 0, 0, 0,
	 BRLOBOL_ACTION_TOGGLE_VIEW_AXES},
	{BRLOBOL_INPUT_KEY_PRESS, '2', BRLOBOL_INPUT_ANY, 0, 0, 0,
	 BRLOBOL_ACTION_VIEW_2},
	{BRLOBOL_INPUT_KEY_PRESS, '3', BRLOBOL_INPUT_ANY, 0, 0, 0,
	 BRLOBOL_ACTION_VIEW_3},
	{BRLOBOL_INPUT_KEY_PRESS, '4', BRLOBOL_INPUT_ANY, 0, 0, 0,
	 BRLOBOL_ACTION_VIEW_4},
	{BRLOBOL_INPUT_KEY_PRESS, '5', BRLOBOL_INPUT_ANY, 0, 0, 0,
	 BRLOBOL_ACTION_VIEW_5},
	{BRLOBOL_INPUT_KEY_PRESS, '6', BRLOBOL_INPUT_ANY, 0, 0, 0,
	 BRLOBOL_ACTION_VIEW_6},
	{BRLOBOL_INPUT_KEY_PRESS, '7', BRLOBOL_INPUT_ANY, 0, 0, 0,
	 BRLOBOL_ACTION_VIEW_7},
	{BRLOBOL_INPUT_KEY_PRESS, 'F', BRLOBOL_INPUT_ANY, 0, 0, 0,
	 BRLOBOL_ACTION_VIEW_FRONT},
	{BRLOBOL_INPUT_KEY_PRESS, 'T', BRLOBOL_INPUT_ANY, 0, 0, 0,
	 BRLOBOL_ACTION_VIEW_TOP},
	{BRLOBOL_INPUT_KEY_PRESS, 'B', BRLOBOL_INPUT_ANY, 0, 0, 0,
	 BRLOBOL_ACTION_VIEW_BOTTOM},
	{BRLOBOL_INPUT_KEY_PRESS, 'L', BRLOBOL_INPUT_ANY, 0, 0, 0,
	 BRLOBOL_ACTION_VIEW_LEFT},
	{BRLOBOL_INPUT_KEY_PRESS, 'R', BRLOBOL_INPUT_ANY,
	 BRLOBOL_INPUT_MOD_SHIFT, 0, 1, BRLOBOL_ACTION_VIEW_REAR},
	{BRLOBOL_INPUT_KEY_PRESS, 'R', BRLOBOL_INPUT_ANY, 0,
	 BRLOBOL_INPUT_MOD_SHIFT, 0, BRLOBOL_ACTION_VIEW_RIGHT},
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

const BRLObolInputProfile &
BRLObolInputContext::defaultViewProfile(void)
{
    return *brlobol_input_default_view_profile();
}
