# smoke -- pipeline test in seconds.  ~30 k cells per quadrant.
# Wall spacing scales too, so y+ is ~200: NOT for results.
H_SCALE = 6.0
# Declared, not inherited: at H_SCALE 6 with H_SCALE_FIN on, ANY tip zone is
# scaled x6 and no longer fits between the tip and ZONE_R[0] -- validate_params
# rejects it and the preset stops working the moment the default changes.  A
# fin tip zone is meaningless at this resolution anyway.
FIN_H_R = None
