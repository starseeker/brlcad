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

/**
 * Toolkit-neutral key values outside the Unicode scalar range.
 *
 * Printable keys continue to use their uppercase Unicode value.  Native
 * adapters translate non-text keys to this range so application action layers
 * do not depend on Qt, Tk, X11, or Win32 key constants.
 */
typedef enum BRLObolInputKey {
    BRLOBOL_INPUT_KEY_F1 = 0x110001,
    BRLOBOL_INPUT_KEY_F2 = 0x110002,
    BRLOBOL_INPUT_KEY_F3 = 0x110003,
    BRLOBOL_INPUT_KEY_F4 = 0x110004,
    BRLOBOL_INPUT_KEY_F5 = 0x110005,
    BRLOBOL_INPUT_KEY_F6 = 0x110006,
    BRLOBOL_INPUT_KEY_F7 = 0x110007,
    BRLOBOL_INPUT_KEY_F8 = 0x110008,
    BRLOBOL_INPUT_KEY_F9 = 0x110009,
    BRLOBOL_INPUT_KEY_F10 = 0x11000a,
    BRLOBOL_INPUT_KEY_F11 = 0x11000b,
    BRLOBOL_INPUT_KEY_F12 = 0x11000c,
    BRLOBOL_INPUT_KEY_SHIFT = 0x110010,
    BRLOBOL_INPUT_KEY_CONTROL = 0x110011,
    BRLOBOL_INPUT_KEY_ALT = 0x110012,
    BRLOBOL_INPUT_KEY_META = 0x110013
} BRLObolInputKey;

/**
 * An action token interpreted by the handler that owns its binding layer.
 *
 * Zero is always BRLOBOL_ACTION_NONE.  The remaining values are deliberately
 * opaque to dispatch: application-owned layers may define any nonzero local
 * identifiers, including values that overlap the standard actions below.
 */
typedef uint32_t BRLObolInputAction;

enum BRLObolStandardInputAction {
    BRLOBOL_ACTION_NONE = 0,
    BRLOBOL_ACTION_VIEW_ROTATE = 1,
    BRLOBOL_ACTION_VIEW_PAN = 2,
    BRLOBOL_ACTION_VIEW_ZOOM = 3,
    BRLOBOL_ACTION_VIEW_CENTER = 4,
    BRLOBOL_ACTION_VIEW_FRONT = 5,
    BRLOBOL_ACTION_VIEW_TOP = 6,
    /* Values 7 and 8 remain unused to preserve existing standard IDs. */
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
};

enum {
    BRLOBOL_INPUT_ANY = -1,
    /* Handler results are intentionally positive when the canvas should
     * schedule a redraw.  Deferred and cancelled results preserve semantic
     * state without making callers infer it from application-local flags. */
    BRLOBOL_INPUT_RESULT_ERROR = -1,
    BRLOBOL_INPUT_RESULT_UNHANDLED = 0,
    BRLOBOL_INPUT_RESULT_HANDLED = 1,
    BRLOBOL_INPUT_RESULT_DEFERRED = 2,
    BRLOBOL_INPUT_RESULT_CANCELLED = 3
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

/**
 * One immutable application-defined action vocabulary and its handler.
 *
 * The context copies the bindings when the layer is installed.  Action IDs
 * are local to this handler and are never classified or interpreted by
 * libbrlobol.
 */
typedef struct BRLOBOL_INPUT_TYPE_EXPORT BRLObolInputActionLayer {
    const char *name;
    const BRLObolInputBinding *bindings;
    size_t bindingCount;
    BRLObolInputActionHandler handler;
} BRLObolInputActionLayer;

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

/**
 * Apply a canonical view-orientation action to a caller-owned bv_context.
 *
 * This handles only the fixed orientation actions shared by Qt, Tk, and
 * MGED.  Editing, selection, snapping, and application-specific faceplate
 * state intentionally remain with their owning application.
 */
BRLOBOL_EXPORT int
brlobol_input_view_orientation_apply(void *view_ctx,
	BRLObolInputAction action);

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

    /**
     * Atomically install one application-owned action layer.  Its bindings
     * win score ties with the base profile and dispatch only to the handler
     * stored in the layer.  An UNHANDLED result falls through to the matching
     * base action.  The owner token guards replacement and teardown.
     */
    int setActionLayer(const BRLObolInputActionLayer *layer, void *owner,
	void *userData);
    int clearActionLayerIf(void *owner);

    int dispatch(const BRLObolInputEvent *event) const;
    int hasAction(BRLObolInputAction action) const;
    size_t bindingCount(void) const;

    static const BRLObolInputProfile &defaultViewProfile(void);

private:
    std::vector<BRLObolInputBinding> bindings;
    std::vector<BRLObolInputBinding> layerBindings;
    BRLObolInputActionHandler handler;
    void *handlerData;
    BRLObolInputActionHandler layerHandler;
    void *layerHandlerData;
    void *layerOwner;
};
#endif

#undef BRLOBOL_INPUT_TYPE_EXPORT

#endif /* BRLOBOL_INPUT_H */
