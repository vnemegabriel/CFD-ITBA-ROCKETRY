#!/usr/bin/env python3
# .venv/bin/python flightModel/post/flowCheck.py <run-dir> [time]   (or pvpython)
#
# Physical checks that do not depend on where the shock is found:
#   entropy     s - s_inf >= 0 in inviscid flow (second law)
#   p0          p0/p0_inf <= 1 there (same law, with h0 constant)
#   h0          h0/h0_inf = 1 there (steady, adiabatic, inviscid)
#   upstream    ahead of the Mach cone of the tip (M > 1) the flow is the free stream
# The first three are tabulated by region and by wall distance, with no cut-off:
# failures that fade with wall distance and with mesh refinement are the
# boundary layer; failures far from the walls, or that do not fall from one
# mesh to the next, are numerical.  The wakes are viscous at any wall distance.
# The fields, the shock sensor Mn = U.grad(p)/(a |grad(p)|) (Mn = 1 on a shock,
# Lovely and Haimes 1999) and the wall distance are written on the plane z = 0
# to <run>/postProcessing/flowCheck/<time>/symmetry.vtp for ParaView.
import sys, math, os
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy as v2n, numpy_to_vtk as n2v

TOL = 1e-3
BANDS = ((0, 0.01), (0.01, 0.03), (0.03, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5), (0.5, 1e3))
GAMMA, R, CP = 1.4, 287.05, 1005.0
WALLS = ('patch/cone', 'patch/fins', 'patch/tail', 'patch/walls')
REGIONS = (('upstream', -1e3, 0.0), ('nose', 0.0, 0.8), ('body', 0.8, 2.529), ('fins', 2.529, 2.830),
           ('boattail', 2.830, 2.955), ('near wake', 2.955, 3.5), ('far wake', 3.5, 1e3))

case = sys.argv[1].rstrip('/')
fc = open(case + '/system/flowConditions').read()
value = lambda k: float(fc.split('\n' + k)[1].split(';')[0])
M, PINF, TINF, ALPHA = value('Minf'), value('pInf'), value('Tinf'), math.radians(value('alpha'))

r = vtk.vtkOpenFOAMReader()
r.SetFileName(case + '/case.foam')
r.UpdateInformation()
times = v2n(r.GetTimeValues())
T = float(sys.argv[2]) if len(sys.argv) > 2 else float(times[-1])
r.SetTimeValue(T)
for i in range(r.GetNumberOfPatchArrays()):
    n = r.GetPatchArrayName(i)
    r.SetPatchArrayStatus(n, int(n == 'internalMesh' or n in WALLS))
r.CreateCellToPointOff()
if hasattr(r, 'DecomposePolyhedraOff'):
    r.DecomposePolyhedraOff()
r.Update()

blocks = {}
it = r.GetOutput().NewIterator()
it.InitTraversal()
while not it.IsDoneWithTraversal():
    blocks[it.GetCurrentMetaData().Get(vtk.vtkCompositeDataSet.NAME())] = it.GetCurrentDataObject()
    it.GoToNextItem()
mesh = blocks['internalMesh']

grad = vtk.vtkGradientFilter()
grad.SetInputData(mesh)
grad.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, 'p')
grad.SetResultArrayName('gradp')
grad.Update()
size = vtk.vtkCellSizeFilter()
size.SetInputData(grad.GetOutput())
size.ComputeVertexCountOff(); size.ComputeLengthOff(); size.ComputeAreaOff(); size.ComputeVolumeOn()
size.Update()
mesh = size.GetOutput()

walls = vtk.vtkAppendPolyData()
for n in WALLS:
    if n.split('/')[1] in blocks:
        g = vtk.vtkGeometryFilter(); g.SetInputData(blocks[n.split('/')[1]]); g.Update()
        walls.AddInputData(g.GetOutput())
walls.Update()
dist = vtk.vtkImplicitPolyDataDistance()
dist.SetInput(walls.GetOutput())
centres = vtk.vtkCellCenters()
centres.SetInputData(mesh)
centres.Update()
sample = vtk.vtkSampleImplicitFunctionFilter()
sample.SetImplicitFunction(dist)
sample.SetInputData(centres.GetOutput())
sample.ComputeGradientsOff()
sample.Update()

cd = mesh.GetCellData()
f = lambda n: v2n(cd.GetArray(n))
x = v2n(centres.GetOutput().GetPoints().GetData())
p, Tc, U, gp = f('p'), f('T'), f('U'), f('gradp')
h = np.cbrt(np.abs(f('Volume')))
wd = np.abs(v2n(sample.GetOutput().GetPointData().GetScalars()))

a = np.sqrt(GAMMA * R * TINF)
UINF = M * a
q = 1 + (GAMMA - 1) / 2 * M**2
P0INF, H0INF = PINF * q**(GAMMA / (GAMMA - 1)), CP * TINF + UINF**2 / 2

ma = np.linalg.norm(U, axis=1) / np.sqrt(GAMMA * R * Tc)
ds = GAMMA / (GAMMA - 1) * np.log(Tc / TINF) - np.log(p / PINF)
p0 = p * (1 + (GAMMA - 1) / 2 * ma**2)**(GAMMA / (GAMMA - 1)) / P0INF
h0 = (CP * Tc + 0.5 * np.einsum('ij,ij->i', U, U)) / H0INF - 1
gn = np.linalg.norm(gp, axis=1)
mn = np.einsum('ij,ij->i', U, gp) / (np.sqrt(GAMMA * R * Tc) * np.maximum(gn, 1e-30))


def where(i):
    return f'x={x[i, 0]:.3f} y={x[i, 1]:.3f} z={x[i, 2]:.3f} m, wall dist {wd[i]:.3f} m, cell {h[i] * 1e3:.1f} mm, Ma {ma[i]:.3f}'


print(f'# {case}  t={T:g}  M={M}  alpha={math.degrees(ALPHA):g} deg  {len(p)} cells  tolerance {TOL}')
heads = ''.join(f'{f"{lo:g}-{hi:g}":>10s}' for lo, hi in BANDS)
for title, v, bad, worst in (
        ('entropy (s-s_inf)/R < -tol', ds, lambda m: ds[m] < -TOL, np.argmin),
        ('p0/p0_inf - 1 > +tol', p0 - 1, lambda m: p0[m] - 1 > TOL, np.argmax),
        ('|h0/h0_inf - 1| > tol', h0, lambda m: np.abs(h0[m]) > TOL, lambda a: np.argmax(np.abs(a)))):
    print(f'\n{title}: % of the cells that fail, by region (rows) and wall distance in m (columns)')
    print(f'  {"region":10s}{heads}')
    for name, xlo, xhi in REGIONS:
        row = ''
        for lo, hi in BANDS:
            m = (x[:, 0] >= xlo) & (x[:, 0] < xhi) & (wd >= lo) & (wd < hi)
            row += f'{100 * np.mean(bad(m)):9.2f}%' if m.sum() >= 20 else f'{"-":>10s}'
        print(f'  {name:10s}{row}')
    print(f'  worst cell in each wall-distance band:')
    for lo, hi in BANDS:
        b = np.where((wd >= lo) & (wd < hi))[0]
        if len(b):
            k = b[worst(v[b])]
            print(f'    {f"{lo:g}-{hi:g}":>10s} m: {v[k]: .3e} at {where(k)}')

if M > 1:
    mu = math.asin(1 / M)
    e = np.array([math.cos(ALPHA), math.sin(ALPHA), 0.0])
    xs = x @ e
    rs = np.linalg.norm(x - np.outer(xs, e), axis=1)
    ahead = (xs < rs / math.tan(mu) - h) & (wd > BANDS[2][0])
    dp = np.abs(p / PINF - 1)
    print(f'\nupstream of the Mach cone of the tip (mu = {math.degrees(mu):.2f} deg): {ahead.sum()} cells')
    for lim in (1e-3, 1e-2):
        b = np.where(ahead & (dp > lim))[0]
        print(f'    |p/p_inf - 1| > {lim:g}: {len(b)} cells', end='')
        if len(b):
            k = b[np.argmax(xs[b] - rs[b] / math.tan(mu))]
            lead = rs[k] / math.tan(mu) - xs[k]
            print(f'; the farthest ahead is {lead:.3f} m ({lead / h[k]:.1f} cells) ahead of the cone at {where(k)}', end='')
        print()

out = os.path.join(case, 'postProcessing', 'flowCheck', f'{T:g}')
os.makedirs(out, exist_ok=True)
for n, v in (('ds_R', ds), ('p0_ratio', p0), ('h0_err', h0), ('Mn', mn), ('wallDist', wd), ('cellSize', h)):
    arr = n2v(np.ascontiguousarray(v, dtype=np.float64), deep=1); arr.SetName(n); cd.AddArray(arr)
plane = vtk.vtkPlane(); plane.SetOrigin(0, 0, 1e-4); plane.SetNormal(0, 0, 1)
cut = vtk.vtkCutter(); cut.SetInputData(mesh); cut.SetCutFunction(plane); cut.Update()
# the appended mode of this VTK writes files its own reader rejects
w = vtk.vtkXMLPolyDataWriter(); w.SetDataModeToBinary()
w.SetFileName(os.path.join(out, 'symmetry.vtp')); w.SetInputData(cut.GetOutput()); w.Write()
print(f'# fields on z = 0: {out}/symmetry.vtp')
