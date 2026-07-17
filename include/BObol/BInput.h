/*                      B I N P U T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BInput.h
 *
 * Toolkit-neutral input normalization and per-owner action dispatch.
 */

#ifndef BOBOL_BINPUT_H
#define BOBOL_BINPUT_H

#include "BObol/BDefines.h"

#include <stdint.h>
#include <stddef.h>
#ifdef __cplusplus
#include <vector>
#  define BOBOL_INPUT_TYPE_EXPORT BOBOL_EXPORT
#else
#  define BOBOL_INPUT_TYPE_EXPORT
#endif

typedef enum BObolInputEventType {
    BOBOL_INPUT_NONE = 0,
    BOBOL_INPUT_KEY_PRESS = 1,
    BOBOL_INPUT_KEY_RELEASE = 2,
    BOBOL_INPUT_POINTER_PRESS = 3,
    BOBOL_INPUT_POINTER_RELEASE = 4,
    BOBOL_INPUT_POINTER_MOTION = 5,
    BOBOL_INPUT_WHEEL = 6,
    BOBOL_INPUT_RESIZE = 7,
    BOBOL_INPUT_FOCUS = 8,
    BOBOL_INPUT_CLOSE = 9,
    BOBOL_INPUT_EXPOSE = 10
} BObolInputEventType;

typedef enum BObolInputModifier {
    BOBOL_INPUT_MOD_NONE = 0u,
    BOBOL_INPUT_MOD_SHIFT = 1u << 0,
    BOBOL_INPUT_MOD_CONTROL = 1u << 1,
    BOBOL_INPUT_MOD_ALT = 1u << 2,
    BOBOL_INPUT_MOD_META = 1u << 3
} BObolInputModifier;

/**
 * Toolkit-neutral key values outside the Unicode scalar range.
 *
 * Printable keys continue to use their uppercase Unicode value.  Native
 * adapters translate non-text keys to this range so application action layers
 * do not depend on Qt, Tk, X11, or Win32 key constants.
 */
typedef enum BObolInputKey {
    BOBOL_INPUT_KEY_F1 = 0x110001,
    BOBOL_INPUT_KEY_F2 = 0x110002,
    BOBOL_INPUT_KEY_F3 = 0x110003,
    BOBOL_INPUT_KEY_F4 = 0x110004,
    BOBOL_INPUT_KEY_F5 = 0x110005,
    BOBOL_INPUT_KEY_F6 = 0x110006,
    BOBOL_INPUT_KEY_F7 = 0x110007,
    BOBOL_INPUT_KEY_F8 = 0x110008,
    BOBOL_INPUT_KEY_F9 = 0x110009,
    BOBOL_INPUT_KEY_F10 = 0x11000a,
    BOBOL_INPUT_KEY_F11 = 0x11000b,
    BOBOL_INPUT_KEY_F12 = 0x11000c,
    BOBOL_INPUT_KEY_SHIFT = 0x110010,
    BOBOL_INPUT_KEY_CONTROL = 0x110011,
    BOBOL_INPUT_KEY_ALT = 0x110012,
    BOBOL_INPUT_KEY_META = 0x110013
} BObolInputKey;

/**
 * An action token interpreted by the handler that owns its binding layer.
 *
 * Zero is always BOBOL_ACTION_NONE.  The remaining values are deliberately
 * opaque to dispatch: application-owned layers may define any nonzero local
 * identifiers, including values that overlap the standard actions below.
 */
typedef uint32_t BObolInputAction;

enum BObolStandardInputAction {
    BOBOL_ACTION_NONE = 0,
    BOBOL_ACTION_VIEW_ROTATE = 1,
    BOBOL_ACTION_VIEW_PAN = 2,
    BOBOL_ACTION_VIEW_ZOOM = 3,
    BOBOL_ACTION_VIEW_CENTER = 4,
    BOBOL_ACTION_VIEW_FRONT = 5,
    BOBOL_ACTION_VIEW_TOP = 6,
    /* Values 7 and 8 remain unused to preserve existing standard IDs. */
    BOBOL_ACTION_VIEW_BOTTOM = 9,
    BOBOL_ACTION_VIEW_LEFT = 10,
    BOBOL_ACTION_VIEW_REAR = 11,
    BOBOL_ACTION_VIEW_RIGHT = 12,
    BOBOL_ACTION_VIEW_2 = 13,
    BOBOL_ACTION_VIEW_3 = 14,
    BOBOL_ACTION_VIEW_4 = 15,
    BOBOL_ACTION_VIEW_5 = 16,
    BOBOL_ACTION_VIEW_6 = 17,
    BOBOL_ACTION_VIEW_7 = 18,
    BOBOL_ACTION_TOGGLE_ADC = 19,
    BOBOL_ACTION_TOGGLE_MODEL_AXES = 20,
    BOBOL_ACTION_TOGGLE_VIEW_AXES = 21,
    BOBOL_ACTION_VIEW_ADJUST = 22,
    BOBOL_ACTION_VIEW_PRIMARY_RELEASE = 23,
    BOBOL_ACTION_VIEW_SECONDARY_RELEASE = 24,
    BOBOL_ACTION_VIEW_CENTER_RELEASE = 25,
    /* Marks the start of a semantic primary-button pan gesture. */
    BOBOL_ACTION_VIEW_PAN_BEGIN = 26,
    /* Marks the end of a semantic primary-button pan gesture. */
    BOBOL_ACTION_VIEW_PAN_END = 27
};

enum {
    BOBOL_INPUT_ANY = -1,
    /* Handler results are intentionally positive when the canvas should
     * schedule a redraw.  Deferred and cancelled results preserve semantic
     * state without making callers infer it from application-local flags. */
    BOBOL_INPUT_RESULT_ERROR = -1,
    BOBOL_INPUT_RESULT_UNHANDLED = 0,
    BOBOL_INPUT_RESULT_HANDLED = 1,
    BOBOL_INPUT_RESULT_DEFERRED = 2,
    BOBOL_INPUT_RESULT_CANCELLED = 3
};

typedef struct BOBOL_INPUT_TYPE_EXPORT BObolInputEvent {
#ifdef __cplusplus
    BObolInputEvent(void);
#endif

    BObolInputEventType type;
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
} BObolInputEvent;

typedef struct BOBOL_INPUT_TYPE_EXPORT BObolInputBinding {
#ifdef __cplusplus
    BObolInputBinding(void);
    BObolInputBinding(BObolInputEventType eventType,
	int key, int button, unsigned int requiredModifiers,
	unsigned int forbiddenModifiers, int priority,
	BObolInputAction action);
#endif

    BObolInputEventType eventType;
    int key;
    int button;
    unsigned int requiredModifiers;
    unsigned int forbiddenModifiers;
    int priority;
    BObolInputAction action;
} BObolInputBinding;

typedef struct BOBOL_INPUT_TYPE_EXPORT BObolInputProfile {
    const char *name;
    const BObolInputBinding *bindings;
    size_t bindingCount;
} BObolInputProfile;

typedef int (*BObolInputActionHandler)(void *userData,
	BObolInputAction action, const BObolInputEvent *event);

/**
 * One immutable application-defined action vocabulary and its handler.
 *
 * The context copies the bindings when the layer is installed.  Action IDs
 * are local to this handler and are never classified or interpreted by
 * libBObol.
 */
typedef struct BOBOL_INPUT_TYPE_EXPORT BObolInputActionLayer {
    const char *name;
    const BObolInputBinding *bindings;
    size_t bindingCount;
    BObolInputActionHandler handler;
} BObolInputActionLayer;

/** Deliver a fully normalized toolkit event to an endpoint-owned context. */
typedef int (*BObolInputEventHandler)(void *userData,
	const BObolInputEvent *event);

__BEGIN_DECLS

/** The standard BRL-CAD view binding profile. */
BOBOL_EXPORT const BObolInputProfile *
bobol_input_default_view_profile(void);

/**
 * Shared keyboard-only view bindings for hosts whose application retains
 * pointer gestures such as edit, selection, or constrained navigation.
 */
BOBOL_EXPORT const BObolInputProfile *
bobol_input_keyboard_view_profile(void);

/**
 * Apply a canonical view-orientation action to a caller-owned bv_context.
 *
 * This handles only the fixed orientation actions shared by Qt, Tk, and
 * MGED.  Editing, selection, snapping, and application-specific faceplate
 * state intentionally remain with their owning application.
 */
BOBOL_EXPORT int
bobol_input_view_orientation_apply(void *view_ctx,
	BObolInputAction action);

__END_DECLS

/**
 * One endpoint- or application-owned binding context.  Neither profiles nor
 * action handlers are process-global, so separate views may safely choose
 * distinct interaction policies.
 */
#ifdef __cplusplus
class BOBOL_EXPORT BObolInputContext {
public:
    BObolInputContext(void);
    ~BObolInputContext(void);

    void clear(void);
    void setProfile(const BObolInputProfile *profile);
    void setBindings(const BObolInputBinding *bindings, size_t count);
    void setActionHandler(BObolInputActionHandler handler, void *userData);
    int clearActionHandlerIf(BObolInputActionHandler handler,
	void *userData);

    /**
     * Atomically install one application-owned action layer.  Its bindings
     * win score ties with the base profile and dispatch only to the handler
     * stored in the layer.  An UNHANDLED result falls through to the matching
     * base action.  The owner token guards replacement and teardown.
     */
    int setActionLayer(const BObolInputActionLayer *layer, void *owner,
	void *userData);
    int clearActionLayerIf(void *owner);

    int dispatch(const BObolInputEvent *event) const;
    int hasAction(BObolInputAction action) const;
    size_t bindingCount(void) const;

    static const BObolInputProfile &defaultViewProfile(void);

private:
    std::vector<BObolInputBinding> bindings;
    std::vector<BObolInputBinding> layerBindings;
    BObolInputActionHandler handler;
    void *handlerData;
    BObolInputActionHandler layerHandler;
    void *layerHandlerData;
    void *layerOwner;
};
#endif

#undef BOBOL_INPUT_TYPE_EXPORT

#endif /* BOBOL_BINPUT_H */
