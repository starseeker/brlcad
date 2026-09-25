/*         D I S P L A Y _ P L A N E _ V I E W _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_DISPLAY_PLANE_VIEW_PRIVATE_H
#define LIBBOBOL_DISPLAY_PLANE_VIEW_PRIVATE_H

#include <Inventor/SbMatrix.h>
#include <Inventor/SbVec2s.h>
#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoCamera.h>

struct BObolDisplayPlaneView {
    SbMatrix worldToClip;
    SbVec2s viewportSize;
};

inline bool
bobol_display_plane_view(const SoViewport &viewport,
	BObolDisplayPlaneView &view)
{
    const SoCamera *camera = viewport.getCamera();
    if (!camera)
	return false;

    view.viewportSize = viewport.getViewportRegion().getViewportSizePixels();
    const float aspect = view.viewportSize[1] > 0 ?
	static_cast<float>(view.viewportSize[0]) / view.viewportSize[1] : 1.0f;
    view.worldToClip = camera->getViewVolume(aspect).getMatrix();
    return true;
}

/** Bind camera projection to one action traversal and always release it. */
class BObolDisplayPlaneViewScope
{
public:
    BObolDisplayPlaneViewScope(const BObolDisplayPlaneView *&slot,
	    const BObolDisplayPlaneView &view) : actionSlot(slot), previous(slot)
    {
	actionSlot = &view;
    }

    ~BObolDisplayPlaneViewScope()
    {
	actionSlot = previous;
    }

    BObolDisplayPlaneViewScope(const BObolDisplayPlaneViewScope &) = delete;
    BObolDisplayPlaneViewScope &operator=(const BObolDisplayPlaneViewScope &) = delete;

private:
    const BObolDisplayPlaneView *&actionSlot;
    const BObolDisplayPlaneView *previous;
};

#endif /* LIBBOBOL_DISPLAY_PLANE_VIEW_PRIVATE_H */
