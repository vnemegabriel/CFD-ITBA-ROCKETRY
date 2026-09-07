# wake_unsteady -- longer, finer near wake for a transient (pimpleFoam) run
# of base-flow shedding.  Near wake to 8 m behind the base, 15 mm cells,
# opened wider so the fine band tracks a spreading, meandering wake.
# Use with --sector full: a symmetry plane suppresses asymmetric shedding.
X_WAKE_2 = 8.00
ZONE0_R_WAKE = 0.90
SEGMENTS = [
    ('up',    dict(h_start=None,   h_end=1.5e-3)),
    ('nose',  dict(h_start=1.5e-3, h_end=0.012)),
    ('cyl',   dict(h_start=0.012,  h_end=0.012)),
    ('tail',  dict(h_start=0.008,  h_end=0.008)),
    ('wake1', dict(h_start=None,   h_end=0.006)),
    ('wake2', dict(h_start=0.006,  h_end=0.015)),
    ('wake3', dict(h_start=0.015,  h_end=1.500)),
]
