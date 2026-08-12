/*              B N A V I G A T I O N G I Z M O . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BNavigationGizmo.h
 *
 * Retained, screen-locked camera orientation control.  The node owns only
 * presentation and renderer-independent hit testing; applications retain
 * view mutation, animation, and input-routing policy.
 */

#ifndef BOBOL_BNAVIGATIONGIZMO_H
#define BOBOL_BNAVIGATIONGIZMO_H

#include "BObol/BDefines.h"

#include <Inventor/SbRotation.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFColor.h>
#include <Inventor/fields/SoSFEnum.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFInt32.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/nodes/SoSeparator.h>

class SoCamera;
class SoFieldSensor;
class SoGLRenderAction;
class SoHUDKit;
class SoSensor;
class SoTranslation;

/**
 * Pixel-space orientation control with cube and axis-ring presentations.
 *
 * Input coordinates passed to hitTest() use the normalized BObol convention:
 * origin at the upper-left, positive X rightward, and positive Y downward.
 * Render geometry uses SoHUDKit's lower-left pixel coordinates internally.
 */
class BOBOL_EXPORT SoBRLNavigationGizmo : public SoSeparator
{
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLNavigationGizmo);

public:
    /** Visual presentation.  Interaction semantics are identical for both. */
    enum Style {
	CUBE = 0,
	CIRCLES = 1
    };

    /** Screen corner occupied by the gizmo. */
    enum Corner {
	LOWER_LEFT = 0,
	LOWER_RIGHT = 1,
	UPPER_LEFT = 2,
	UPPER_RIGHT = 3
    };

    /**
     * Pickable orientation part.  Direction names encode outward world-space
     * directions; cube faces have one component, edges two, and corners
     * three.  PART_ORBIT is the circles style's non-snapping drag field.
     */
    enum Part {
	PART_NONE = 0,
	PART_NEG_X_NEG_Y_NEG_Z = 1,
	PART_NEG_X_NEG_Y = 2,
	PART_NEG_X_NEG_Y_POS_Z = 3,
	PART_NEG_X_NEG_Z = 4,
	PART_NEG_X = 5,
	PART_NEG_X_POS_Z = 6,
	PART_NEG_X_POS_Y_NEG_Z = 7,
	PART_NEG_X_POS_Y = 8,
	PART_NEG_X_POS_Y_POS_Z = 9,
	PART_NEG_Y_NEG_Z = 10,
	PART_NEG_Y = 11,
	PART_NEG_Y_POS_Z = 12,
	PART_NEG_Z = 13,
	PART_POS_Z = 15,
	PART_POS_Y_NEG_Z = 16,
	PART_POS_Y = 17,
	PART_POS_Y_POS_Z = 18,
	PART_POS_X_NEG_Y_NEG_Z = 19,
	PART_POS_X_NEG_Y = 20,
	PART_POS_X_NEG_Y_POS_Z = 21,
	PART_POS_X_NEG_Z = 22,
	PART_POS_X = 23,
	PART_POS_X_POS_Z = 24,
	PART_POS_X_POS_Y_NEG_Z = 25,
	PART_POS_X_POS_Y = 26,
	PART_POS_X_POS_Y_POS_Z = 27,
	/** Draggable ring field with no canonical click orientation. */
	PART_ORBIT = 28
    };

    SoSFString overlayId; /**< Stable feature identity. */
    SoSFBool visible; /**< Master presentation visibility. */
    SoSFEnum style; /**< SoBRLNavigationGizmo::Style. */
    SoSFEnum corner; /**< SoBRLNavigationGizmo::Corner. */
    SoSFFloat size; /**< Overall control size in physical pixels. */
    SoSFFloat margin; /**< Distance from the selected viewport corner. */
    SoSFColor xColor; /**< X-axis color. */
    SoSFColor yColor; /**< Y-axis color. */
    SoSFColor zColor; /**< Z-axis color. */
    SoSFColor panelColor; /**< Translucent panel color. */
    SoSFColor highlightColor; /**< Hover/active accent color. */
    SoSFInt32 hoverPart; /**< Current hover Part. */
    SoSFInt32 activePart; /**< Current pressed Part. */

    SoBRLNavigationGizmo(void);
    static void initClass(void);

    /** Attach a borrowed camera and track its orientation until replaced. */
    void setCamera(SoCamera *camera);

    /** Return the borrowed tracked camera, or NULL. */
    SoCamera *getCamera(void) const;

    /** Rebuild the small retained HUD subtree from the current fields. */
    SoHUDKit *rebuildGeometry(void);

    /** Return the current HUD kit, or NULL when not realized. */
    SoHUDKit *getHUDKit(void) const;

    /** Update hover feedback and rebuild only when the part changed. */
    void setHoverPart(Part part);

    /** Update pressed feedback and rebuild only when the part changed. */
    void setActivePart(Part part);

    /**
     * Resolve a normalized upper-left-origin pixel to an interactive part.
     * Width and height are the physical-pixel viewport dimensions.
     */
    Part hitTest(int x, int y, int width, int height) const;

    /** Decode a part to its {-1,0,+1} world direction. */
    static SbBool partDirection(Part part, SbVec3f &direction);

    /** Convert a part direction to BRL-CAD azimuth/elevation degrees. */
    static SbBool partAet(Part part, float &azimuth, float &elevation);

    /** Keep the screen anchor synchronized with the active render viewport. */
    virtual void GLRender(SoGLRenderAction *action) override;

protected:
    virtual ~SoBRLNavigationGizmo(void);

private:
    static void cameraChanged(void *userData, SoSensor *sensor);
    void syncCamera(void);
    void updateAnchor(int width, int height);
    SbRotation displayRotation(void) const;

    SoCamera *trackedCamera;
    SoFieldSensor *cameraSensor;
    SoHUDKit *hudKit;
    SoTranslation *anchorTranslation;
    SbRotation rotation;
};

#endif /* BOBOL_BNAVIGATIONGIZMO_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
