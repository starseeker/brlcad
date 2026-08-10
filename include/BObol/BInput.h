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

/** @addtogroup bobol_query
 * @{ */

/** @ingroup bobol_query
 * Normalized event kind.  Values are toolkit-independent and may be passed
 * between a host adapter and its owner-thread endpoint dispatcher.
 */
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

/** @ingroup bobol_query
 * Bit mask of normalized keyboard modifiers.
 */
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

/** @ingroup bobol_query
 * Standard BRL-CAD view actions understood by the default profile.
 * Application action layers may use other nonzero identifiers.
 */
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

/** @ingroup bobol_query
 * One borrowed, immutable event delivered synchronously on the endpoint owner
 * thread.  Pixel coordinates use the bound viewport with its native origin;
 * deltas and resize dimensions are pixels, and timestamp units are supplied by
 * the host adapter.
 */
typedef struct BOBOL_INPUT_TYPE_EXPORT BObolInputEvent {
#ifdef __cplusplus
    BObolInputEvent(void);
#endif

    BObolInputEventType type; /**< Event kind. */
    uint64_t timestamp; /**< Host-supplied monotonic timestamp. */
    int x; /**< Current viewport x coordinate in pixels. */
    int y; /**< Current viewport y coordinate in pixels. */
    int dx; /**< Horizontal motion delta in pixels. */
    int dy; /**< Vertical motion delta in pixels. */
    int wheelDelta; /**< Signed wheel step or high-resolution delta. */
    int button; /**< Changed button number, or zero when not applicable. */
    unsigned int buttons; /**< Host-defined mask of currently pressed buttons. */
    int key; /**< Unicode scalar or BObolInputKey value. */
    unsigned int modifiers; /**< BObolInputModifier bit mask. */
    unsigned int width; /**< Resize width in physical pixels. */
    unsigned int height; /**< Resize height in physical pixels. */
} BObolInputEvent;

/** @ingroup bobol_query
 * Immutable event-to-action rule copied into a BObolInputContext.  A field set
 * to BOBOL_INPUT_ANY is a wildcard; higher priority wins before specificity
 * and installation order are considered.
 */
typedef struct BOBOL_INPUT_TYPE_EXPORT BObolInputBinding {
#ifdef __cplusplus
    BObolInputBinding(void);
    BObolInputBinding(BObolInputEventType eventType,
	int key, int button, unsigned int requiredModifiers,
	unsigned int forbiddenModifiers, int priority,
	BObolInputAction action);
#endif

    BObolInputEventType eventType; /**< Event kind to match. */
    int key; /**< Key value to match, or BOBOL_INPUT_ANY. */
    int button; /**< Button to match, or BOBOL_INPUT_ANY. */
    unsigned int requiredModifiers; /**< Modifier bits that must be present. */
    unsigned int forbiddenModifiers; /**< Modifier bits that must be absent. */
    int priority; /**< Tie-breaking priority; larger values win. */
    BObolInputAction action; /**< Opaque action delivered to the owning handler. */
} BObolInputBinding;

/** @ingroup bobol_query
 * Borrowed immutable binding profile.  Setters copy its bindings; the name
 * and arrays need only remain valid for the duration of the setter call.
 */
typedef struct BOBOL_INPUT_TYPE_EXPORT BObolInputProfile {
    const char *name; /**< Diagnostic profile name. */
    const BObolInputBinding *bindings; /**< Binding array, or NULL when empty. */
    size_t bindingCount; /**< Number of entries in bindings. */
} BObolInputProfile;

/** @ingroup bobol_query
 * Synchronous owner-thread action callback.  @p event is borrowed for the
 * call.  Return a BOBOL_INPUT_RESULT_* value; callbacks may queue work but
 * must not destroy or recursively dispatch through their endpoint.
 */
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
    const char *name; /**< Diagnostic layer name, borrowed during installation. */
    const BObolInputBinding *bindings; /**< Rules copied during installation. */
    size_t bindingCount; /**< Number of rules in bindings. */
    BObolInputActionHandler handler; /**< Required handler for this vocabulary. */
} BObolInputActionLayer;

/** @ingroup bobol_query
 * Deliver a fully normalized toolkit event to an endpoint-owned context.
 * Both pointers are borrowed for the synchronous owner-thread call; return a
 * BOBOL_INPUT_RESULT_* value.
 */
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

/** @} */

#endif /* BOBOL_BINPUT_H */
