# fintip -- spanwise refinement at the fin tip (r = 0.2355 m).
# FIN_H_R inserts a zone boundary at the outer edge of the tip smear with a
# 6 mm radial cell, so the fin closes on a node line and the tip vortex is
# resolved: 6 mm at the tip instead of 12.  Combine with any other preset:
#     python build.py --preset medium --preset fintip
FIN_H_R = 0.006
