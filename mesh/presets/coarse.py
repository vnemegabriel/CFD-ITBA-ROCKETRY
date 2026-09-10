# coarse -- first look at a flow field.  ~473 k cells per quadrant, ~45 s.
# Keeps y1 (y+ = 32) at the wall and the FIN at full resolution: 5 mm along
# the chord, 0.5 mm normal to it, 12 mm spanwise at the tip -- which is also
# the height of the tip chamfer, since that is max(FIN_TIP_SMEAR, the radial
# cell there).  `--preset coarse --preset fintip` takes it to 6 mm.  Everything
# else is 3x coarser than `fine`.  ZONE_H[0] is lowered so the step from the
# 12 mm fin zone to zone 0 stays under 1.3 per cell.
H_SCALE = 3.0
H_SCALE_WALL = False
H_SCALE_FIN = False
FIN_H_R = 0.012
ZONE_H = [0.010, 0.200, 1.600]
