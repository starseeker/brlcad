/*          W I N D O W _ H O S T _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_WINDOW_HOST_PRIVATE_H
#define LIBBOBOL_WINDOW_HOST_PRIVATE_H

#include "BObol/BViewController.h"
#include "BObol/BWindowHost.h"
#include "view_controller_private.h"

#include <memory>

class SoGroup;

struct BObolWindowHost::PreparedOpenPublication {
    std::unique_ptr<BObolWindowDesc> desc;
    BObolViewController::PreparedViewportPublication viewport;
    bool changed = false;
};

SoGroup *bobol_window_host_root_group(BObolViewController *controller);
void bobol_window_host_detach_controller_noexcept(
    BObolWindowHost *host) noexcept;

#endif /* LIBBOBOL_WINDOW_HOST_PRIVATE_H */
