# fintip -- spanwise refinement at the fin tip (r = 0.2355 m).
# FIN_H_R inserts a zone boundary AT FIN_TIP_R with a 6 mm radial cell, so the
# fin closes on a node line and the tip vortex is resolved: 6 mm at the tip
# instead of 12.  It halves the tip chamfer with it, because the chamfer is
# max(FIN_TIP_SMEAR, the radial cell at the tip).  Combine with any other preset:
#     python build.py --preset medium --preset fintip
FIN_H_R = 0.006
