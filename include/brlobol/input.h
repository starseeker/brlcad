/*                         I N P U T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/input.h
 *
 * Toolkit-neutral input normalization and per-owner action dispatch.
 */

#ifndef BRLOBOL_INPUT_H
#define BRLOBOL_INPUT_H

#include "brlobol/defines.h"

#include <stdint.h>
#include <stddef.h>
#ifdef __cplusplus
#include <vector>
#  define BRLOBOL_INPUT_TYPE_EXPORT BRLOBOL_EXPORT
#else
#  define BRLOBOL_INPUT_TYPE_EXPORT
#endif

typedef enum BRLObolInputEventType {
    BRLOBOL_INPUT_NONE = 0,
    BRLOBOL_INPUT_KEY_PRESS = 1,
    BRLOBOL_INPUT_KEY_RELEASE = 2,
    BRLOBOL_INPUT_POINTER_PRESS = 3,
    BRLOBOL_INPUT_POINTER_RELEASE = 4,
    BRLOBOL_INPUT_POINTER_MOTION = 5,
    BRLOBOL_INPUT_WHEEL = 6,
    BRLOBOL_INPUT_RESIZE = 7,
    BRLOBOL_INPUT_FOCUS = 8,
    BRLOBOL_INPUT_CLOSE = 9,
    BRLOBOL_INPUT_EXPOSE = 10
} BRLObolInputEventType;

typedef enum BRLObolInputModifier {
    BRLOBOL_INPUT_MOD_NONE = 0u,
    BRLOBOL_INPUT_MOD_SHIFT = 1u << 0,
    BRLOBOL_INPUT_MOD_CONTROL = 1u << 1,
    BRLOBOL_INPUT_MOD_ALT = 1u << 2,
    BRLOBOL_INPUT_MOD_META = 1u << 3
} BRLObolInputModifier;

typedef enum BRLObolInputAction {
    BRLOBOL_ACTION_NONE = 0,
    BRLOBOL_ACTION_VIEW_ROTATE = 1,
    BRLOBOL_ACTION_VIEW_PAN = 2,
    BRLOBOL_ACTION_VIEW_ZOOM = 3,
    BRLOBOL_ACTION_VIEW_CENTER = 4,
    BRLOBOL_ACTION_VIEW_FRONT = 5,
    BRLOBOL_ACTION_VIEW_TOP = 6,
    BRLOBOL_ACTION_APP_CANCEL = 7,
    BRLOBOL_ACTION_APP_ACCEPT = 8,
    BRLOBOL_ACTION_VIEW_BOTTOM = 9,
    BRLOBOL_ACTION_VIEW_LEFT = 10,
    BRLOBOL_ACTION_VIEW_REAR = 11,
    BRLOBOL_ACTION_VIEW_RIGHT = 12,
    BRLOBOL_ACTION_VIEW_2 = 13,
    BRLOBOL_ACTION_VIEW_3 = 14,
    BRLOBOL_ACTION_VIEW_4 = 15,
    BRLOBOL_ACTION_VIEW_5 = 16,
    BRLOBOL_ACTION_VIEW_6 = 17,
    BRLOBOL_ACTION_VIEW_7 = 18,
    BRLOBOL_ACTION_TOGGLE_ADC = 19,
    BRLOBOL_ACTION_TOGGLE_MODEL_AXES = 20,
    BRLOBOL_ACTION_TOGGLE_VIEW_AXES = 21,
    BRLOBOL_ACTION_VIEW_ADJUST = 22,
    BRLOBOL_ACTION_VIEW_PRIMARY_RELEASE = 23,
    BRLOBOL_ACTION_VIEW_SECONDARY_RELEASE = 24,
    BRLOBOL_ACTION_VIEW_CENTER_RELEASE = 25,
    /* Marks the start of a semantic primary-button pan gesture. */
    BRLOBOL_ACTION_VIEW_PAN_BEGIN = 26,
    /* Marks the end of a semantic primary-button pan gesture. */
    BRLOBOL_ACTION_VIEW_PAN_END = 27
} BRLObolInputAction;

enum {
    BRLOBOL_INPUT_ANY = -1
};

typedef struct BRLOBOL_INPUT_TYPE_EXPORT BRLObolInputEvent {
#ifdef __cplusplus
    BRLObolInputEvent(void);
#endif

    BRLObolInputEventType type;
    uint64_t timestamp;
    int x;
    int y;
    int dx;
    int dy;
    int wheelDelta;
    int button;
    unsigned int buttons;
    int key;
    unsigned int modifiers;
    unsigned int width;
    unsigned int height;
} BRLObolInputEvent;

typedef struct BRLOBOL_INPUT_TYPE_EXPORT BRLObolInputBinding {
#ifdef __cplusplus
    BRLObolInputBinding(void);
    BRLObolInputBinding(BRLObolInputEventType eventType,
	int key, int button, unsigned int requiredModifiers,
	unsigned int forbiddenModifiers, int priority,
	BRLObolInputAction action);
#endif

    BRLObolInputEventType eventType;
    int key;
    int button;
    unsigned int requiredModifiers;
    unsigned int forbiddenModifiers;
    int priority;
    BRLObolInputAction action;
} BRLObolInputBinding;

typedef struct BRLOBOL_INPUT_TYPE_EXPORT BRLObolInputProfile {
    const char *name;
    const BRLObolInputBinding *bindings;
    size_t bindingCount;
} BRLObolInputProfile;

typedef int (*BRLObolInputActionHandler)(void *userData,
	BRLObolInputAction action, const BRLObolInputEvent *event);

/** Deliver a fully normalized toolkit event to an endpoint-owned context. */
typedef int (*BRLObolInputEventHandler)(void *userData,
	const BRLObolInputEvent *event);

__BEGIN_DECLS

/** The standard BRL-CAD view binding profile. */
BRLOBOL_EXPORT const BRLObolInputProfile *
brlobol_input_default_view_profile(void);

/**
 * Shared keyboard-only view bindings for hosts whose application retains
 * pointer gestures such as edit, selection, or constrained navigation.
 */
BRLOBOL_EXPORT const BRLObolInputProfile *
brlobol_input_keyboard_view_profile(void);

__END_DECLS

/**
 * One endpoint- or application-owned binding context.  Neither profiles nor
 * action handlers are process-global, so separate views may safely choose
 * distinct interaction policies.
 */
#ifdef __cplusplus
class BRLOBOL_EXPORT BRLObolInputContext {
public:
    BRLObolInputContext(void);
    ~BRLObolInputContext(void);

    void clear(void);
    void setProfile(const BRLObolInputProfile *profile);
    void setBindings(const BRLObolInputBinding *bindings, size_t count);
    void setActionHandler(BRLObolInputActionHandler handler, void *userData);
    int clearActionHandlerIf(BRLObolInputActionHandler handler,
	void *userData);

    int dispatch(const BRLObolInputEvent *event) const;
    int hasAction(BRLObolInputAction action) const;
    size_t bindingCount(void) const;

    static const BRLObolInputProfile &defaultViewProfile(void);

private:
    std::vector<BRLObolInputBinding> bindings;
    BRLObolInputActionHandler handler;
    void *handlerData;
};
#endif

#undef BRLOBOL_INPUT_TYPE_EXPORT

#endif /* BRLOBOL_INPUT_H */
