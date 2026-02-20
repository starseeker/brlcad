# setup_sim.tcl - Create the truck-on-terrain simulation database
# 
# Run from share/db directory (so terra.bin is accessible):
#   cd /path/to/brlcad-build/share/db
#   mged -c output.g < setup_sim.tcl
#
# Geometry facts (all in mm):
#   Truck (m35 component) bounding box:
#     min: {-1827, -1250, 0}  max: {4892, 1250, 2618}
#     center X ~ 1532, center Y = 0
#   Terrain center: X=12800000, Y=12800000
#   Terrain center elevation: ~1237625 mm
#   Truck placed 500m (500000 mm) above terrain center elevation
#
# Translation to apply to truck:
#   dx = 12800000 - 1532 = 12798468
#   dy = 12800000 - 0    = 12800000
#   dz = 1237625 + 500000 = 1737625

dbconcat terra.g terra_
dbconcat m35.g m35_

# Static terrain body (mass=0 means it won't move)
comb terrain.sim u terra_ground.r
attr set terrain.sim simulate::type region
attr set terrain.sim simulate::mass 0.0
attr set terrain.sim simulate::roi_proxy 1

# Dynamic truck body (2000 kg)
comb truck.sim u m35_component
attr set truck.sim simulate::type region
attr set truck.sim simulate::mass 2000.0

# Scene combination with gravity
comb scene.c u terrain.sim u truck.sim
attr set scene.c simulate::gravity {<0,0,-9.80665>}

# Position truck above terrain center using matrix on truck.sim in scene.c
# The matrix is: identity rotation, with translation (dx, dy, dz)
arced scene.c/truck.sim matrix rarc 1 0 0 12798468 0 1 0 12800000 0 0 1 1737625 0 0 0 1
