# Mesh

This document gives the mesh procedure, the requirements for the geometry,
and the dictionaries in `mesh/`. [WORKFLOW.md](WORKFLOW.md) §3 gives the
command.

## 1. Mesh generator

The mesh generator is cfMesh (`cartesianMesh`). It makes an octree mesh from
one maximum cell size. All other sizes are refinement levels. Each level
divides the cell size by 2:

```
cell size = maxCellSize / 2^level
```

With `maxCellSize 1.6`:

| Level | Cell size |
|---|---|
| 2 | 400 mm |
| 4 | 100 mm |
| 5 | 50 mm |
| 6 | 25 mm |
| 7 | 12.5 mm |
| 8 | 6.25 mm |
| 9 | 3.13 mm |
| 10 | 1.56 mm |

A level has an effect only on its patch or its region. If a level is
incorrect, the cells in that region are too large or too small. No other
region changes.

## 2. Procedure in Allmesh

`Allmesh` does these steps in the mesh directory. If an application writes
`FOAM FATAL` in its log, `Allmesh` shows the error and stops.

| Step | Application | Input | Output | Function |
|---|---|---|---|---|
| 1 | `cp` | `mesh/system/`, `mesh/meshInfo` | `<mesh-dir>/system/`, `<mesh-dir>/meshInfo` | Each mesh keeps a copy of its dictionaries. |
| 2 | `sed`, `foamDictionary` | `meshInfo`, `--maxCellSize`, `--nLayers`, `--thicknessRatio` | `system/meshDict` | Writes `lBody` and `rBody` of `meshInfo` (§4.3), then the options, into the copy. |
| 2a | `cp` | `<stl>_finEdges.obj`, `<stl>_finRoots.obj` | `finEdges.obj`, `finRoots.obj` | The fin edges and fin roots for `edgeMeshRefinement`. Refer to §4.5. |
| 3 | `sed` | STL | `surface.stl` | Changes the solid names to the patch names. Refer to §3. |
| 4 | `surfaceGenerateBoundingBox` | `surface.stl` | `geometry.stl` | Adds the six faces of the domain, symmetric about z = 0. Refer to §5. |
| 5 | `surfaceFeatureEdges -angle 60` | `geometry.stl` | `geometry.fms` | Finds the feature edges. Refer to §2.1. |
| 6 | `cartesianMesh` | `geometry.fms`, `system/meshDict` | `constant/polyMesh` | Makes the mesh. |
| 7 | `createPatch -overwrite` | `createPatchDict` | `constant/polyMesh` | Merges the six domain faces into `inlet`, `outlet` and `box`. |
| 7a | `topoSet`, `subsetMesh upper -overwrite` | `topoSetDictHalf` | `constant/polyMesh` | Half mesh only: keeps the cells with the centre at z ≥ 0. The cut faces go to `oldInternalFaces`. Refer to §5.1. |
| 7b | `mesh/snapPlane.py . oldInternalFaces` | `constant/polyMesh/points` | `constant/polyMesh/points` | Half mesh only: sets z = 0 on the cut. |
| 7c | `createPatch -overwrite` | `createPatchDictHalf` | `constant/polyMesh` | Half mesh only: `oldInternalFaces` becomes `symm`, of type `symmetry`. Sets `sector`, `symmetryFraction` and `Aref` in `meshInfo`. |
| 8 | `checkMesh -allTopology -allGeometry` | mesh | `log.checkMesh` | Measures the mesh quality. |
| 9 | `foamDictionary` | `log.checkMesh`, `meshDict` | `meshInfo` | Writes `nCells`, `maxCellSize`, `finLevel`, `nLayers`, `thicknessRatio`, `meshDir`. |
| 10 | `grep` | `log.checkMesh` | screen | Shows the summary. Stops if there are negative cell volumes. |

### 2.1 Function of the .fms file

An STL file contains only triangles and solid names. It does not contain
edges. The `.fms` file contains the triangles, the patches and the feature
edges.

`surfaceFeatureEdges` marks an edge as a feature edge when the angle between
the two adjacent triangles is more than 60°. `cartesianMesh` then moves the
mesh vertices onto these edges.

The rocket has no sharp edge left: the fin edges, the fin root, the shoulder
and the base edge are rounded (§3.3), and each round is made of facets of
less than 30°. With 30°, small triangles where the root fillet meets the
body gave 172 false feature edges (56 at the foot of the leading-edge
fillet, 64 at the fin tips), and `cartesianMesh` bent the cells along them.
With 60° only the edges of the domain box (12) and the nose (88) remain.

`cartesianMesh` also accepts an STL file in `surfaceFile`. In that case, the
mesh has no feature edges.

## 3. Geometry requirements

| Item | Requirement | Reason |
|---|---|---|
| Format | ASCII STL | A binary STL has no solid names. All the rocket then becomes one patch. |
| Solid names | `nosecone`, `body`, `boattail`, `fins` | `Allmesh` changes the first three to `cone`, `walls`, `tail`. The dictionaries and the case use these four names. |
| Units | Metres | `meshDict` and `meshInfo` use metres. |
| Position | Nose tip at the origin, axis along +x | The refinement regions in `meshDict`, the fin region in `meshInfo` and the moment reference `CofR (0 0 0)` use this position. |

### 3.1 Solid names that are not in the list

`Allmesh` does not stop if a solid has a different name. The patch keeps the
name of the solid. These are the effects:

| Step | Effect |
|---|---|
| `cartesianMesh`, `localRefinement` | An entry without a patch gives only a warning. The patch gets no surface refinement. The cell size near the wall is then the size of the volume region (50 mm in `nearField`), not 6.25 mm. |
| `cartesianMesh`, `boundaryLayers` | The patch gets the global `nLayers 10`. On a thin fin, 10 layers can cause tangled cells. |
| Case, `0.orig/` | The boundary conditions use the names `cone`, `walls`, `tail`, `fins`. A patch with a different name has no boundary condition, and the solver stops. |
| Case, `forceCoeffs` | `wallPatches` in `meshInfo` does not contain the patch. Its force is not in the coefficients. |

Before you make a mesh from a new STL, examine the solid names:

```bash
grep -E '^ *solid ' mesh/stl/<file>.stl
```

### 3.2 Use a new geometry

1. Put the STL in `mesh/stl/`, or write it with `makeRocket.py` (§3.3).
   The `_finEdges.obj` file must be next to it.
2. Make sure that the STL agrees with the requirements in §3.
3. Change the distances in `surfaceGenerateBoundingBox` in `Allmesh`. Refer
   to §5.
4. Write `lBody` and `rBody` of the new rocket in `mesh/meshInfo`. The
   regions in `objectRefinements` scale with them. Refer to §4.3.
5. Change `mesh/meshInfo`. Refer to §7.
6. Run `mesh/Allmesh <mesh-dir> mesh/stl/<file>.stl`.

### 3.3 makeRocket.py

`mesh/makeRocket.py` writes the geometry from its dimensions. Run it with no
option for the Aconcagua:

```bash
python3 mesh/makeRocket.py mesh/stl/Aconcagua_biconvex117.stl
```

It writes two files: the STL and `<name>_finEdges.obj`, the leading and
trailing edges of the fins as lines. `python3 mesh/makeRocket.py -h` gives
the options. A file with one `name value` on each line can give the
dimensions (`--dims`).

The geometry of the half meshes (§5.1) comes from the same script:

```bash
python3 mesh/makeRocket.py mesh/stl/Aconcagua_biconvex117_roll45.stl --finRoll 45
```

`--finRoll 45` turns the fins to 45, 135, 225 and 315°, so that the plane
z = 0 passes between two fins. The STL is the whole rocket: the half mesh is
cut from the full one (§5.1).

| Part | Shape | Solid name |
|---|---|---|
| Nose | Von Kármán, 0.8 m long | `nosecone` |
| Body | Cylinder, d = 0.151 m | `body` |
| Fins | 4 fins, biconvex section of constant thickness 11.7 mm, edges rounded with 1 mm, a 3 mm fillet all round the root | `fins` |
| Boattail | Cone to d = 0.110 m, and the base; shoulder and base edge rounded with 3 mm | `boattail` |

The script writes a file only if the surface is closed and all the normals
point out. These properties of the surface come from measurements with
`cartesianMesh`. Do not change them:

| Property | Reason, measured with the same `meshDict` |
|---|---|
| Nose and body + fins + boattail are two closed shells that touch on coplanar caps. The rings of the two shells are turned by half a facet. | With a sharp shoulder, one continuous hull gave 12 to 87 negative cell volumes at the body-boattail shoulder, so the boattail was a third shell. With the shoulder rounded (`--shoulderRadius`, 3 mm) the continuous hull gave 0 on `h1`. The boattail must be in the body shell, because the root fillet of the trailing edge runs 3 mm onto the cone. |
| Body and fins are one shell: the fin root is cut into the cylinder. | With the fins as separate shells sunk into the body, the fin-body junction was not an edge of the surface. The layers of the two surfaces met at one node in the corner. |
| The leading edge, the trailing edge and the tip of the fins have a radius (`--finEdgeRadius`, 1 mm). | With sharp edges, the layers of the two fin faces met at the edge: maximum skewness 25 to 33, non-orthogonality up to 176°. With the radius: skewness 3.7, non-orthogonality 77° to 82°. |
| A fillet of 3 mm all round each fin root (`--finRootRadius`), from the fin surface at a height of 3 mm to the body 3 mm off the root. | Without it the fin face and the body met at about 90°, and the 3 fin layers and the 10 body layers closed to one node along the root. With it, on `h1`: non-orthogonality 81.9° → 67.7° (0 faces above 70°, 156 before), skewness 3.58 → 1.89, faces with interpolation weight below 0.05 478 → 250, smallest cell 2.6 → 6.9 × 10⁻¹³ m³, 1 577 819 → 1 583 757 cells. Refinement spheres at the root corners did not help: the corner was already at the level of `finEdges`. |
| The base edge is rounded with 3 mm (`--baseEdgeRadius`). | |
| The root trailing edge is on the shoulder (`--finFromBase` = `--tailLength`). | |
| 44 facets around the axis (`--nTheta`). | With 96 facets the continuous hull had more negative volumes. |

The radius is in the section parallel to the axis. Normal to the swept
leading edge (57°), it is approximately 0.5 mm. The section is the convex
hull of the biconvex arcs and a circle of the radius at each edge, so the
profile is a little thicker in approximately the first 80 mm of the chord.

## 4. meshDict

File: `mesh/system/meshDict`. Each mesh directory has its own copy.

### 4.1 Global entries

| Entry | Value | Function |
|---|---|---|
| `surfaceFile` | `"geometry.fms"` | Surface that `Allmesh` makes in step 5. |
| `maxCellSize` | `1.6` | Maximum cell size in m. All levels refer to this size. `Allmesh --maxCellSize` changes it. |

### 4.2 localRefinement

Surface refinement for each patch. `refinementThickness` is the distance
from the wall in which the cells have the level of the patch. Outside this
distance, the cell size increases to `maxCellSize` one level at a time.

| Patch | Level | Cell size | `refinementThickness` | Reason |
|---|---|---|---|---|
| `cone` | 8 | 6.25 mm | 0.060 m | The bow shock is in this band. |
| `walls` | 8 | 6.25 mm | 0.050 m | Body boundary layer. |
| `tail` | 8 | 6.25 mm | 0.050 m | The expansion on the boattail sets the base pressure. |
| `fins` | 9 | 3.13 mm | 0.020 m | Fin chord. |

### 4.3 objectRefinements

Volume refinement in regions. A `box` has a centre and three lengths. A
`cone` has two points on the axis and a radius at each point.

| Region | Type | Level | Cell size | Extent |
|---|---|---|---|---|
| `shoulder` | hollowCone | 9 | 3.13 mm | x from 2.81 to 2.85 m, r from 0.065 to 0.090 m: the body-boattail shoulder |
| `baseEdge` | hollowCone | 10 | 1.56 mm | x from 2.94 to 2.97 m, r from 0.045 to 0.065 m: the base edge, where the base shear layer starts |
| `nearField` | box | 5 | 50 mm | x from −0.20 to 3.20 m, y and z ±0.50 m |
| `baseWake` | cone | 8 | 6.25 mm | x from 2.94 to 3.30 m, radius 0.08 m |
| `nearWake` | cone | 7 | 12.5 mm | x from 2.90 to 3.80 m, radius 0.14 m |
| `midWake` | cone | 6 | 25 mm | x from 2.85 to 5.00 m, radius 0.26 m |
| `wake` | box | 4 | 100 mm | x from 1.0 to 9.0 m, y and z ±0.80 m |
| `farWake` | box | 2 | 400 mm | x from 1.0 to 27.0 m, y and z ±2.0 m |
| `noseTip` | sphere | 9 | 3.13 mm | Centre (0 0 0), radius 0.040 m |
| `noseApex` | sphere | 11 | 0.78 mm | Centre (0 0 0), radius 0.010 m |

`noseTip` has the level of the fins on each mesh. `noseApex` makes the curve
of the tip. The layer thickness follows the surface cell size, so it also
keeps the layers thinner than the radius of the nose near the tip. Without
`noseApex`, the full mesh `h2` with t₁ = 0.241 mm had 13 negative cell
volumes at the tip.

On the full `h1`, `noseApex` gave cells of 4 µm at the apex. With
`rhoCentralFoam` the smallest cell sets the time step of the whole mesh.
`noseApex` and `finEdges` at the fin level, with t₁ = 0.7 mm, raised it to
8.8 µm, but the surface of the fins lost a third of its faces. The half mesh
keeps the levels of the full one: it is the full mesh cut in two.

The wake regions decrease the cell size one level at a time behind the base.
The numbers are for the Aconcagua (length 2.955 m, body radius 0.0755 m).
Each region is written as `#eval{ <value>*$sx }` along x and
`#eval{ <value>*$sr }` across, with `sx = lBody/2.955` and
`sr = rBody/0.0755`. `Allmesh` writes `lBody` and `rBody` from `meshInfo`
into the copy of `meshDict`, so a rocket of another size gets the regions
moved and scaled with it, and the Aconcagua gets the values of the table.

> **CAUTION** `foamDictionary -set` evaluates the `#eval` entries and writes
> the numbers. `Allmesh` writes `lBody` and `rBody` with `sed` before the
> first `foamDictionary -set`, so the numbers are those of `meshInfo`. Do not
> set `lBody` with `foamDictionary` in a copy that was already written.

### 4.4 boundaryLayers

cfMesh extrudes the layers from the surface mesh. The first cell thickness
thus comes from the surface cell size, `nLayers` and `thicknessRatio`. The
dictionary does not set a y+ target.

| Entry | Value | Function |
|---|---|---|
| `nLayers` | 1 | Layers on all patches that have no entry in `patchBoundaryLayers`: the six faces of the domain. One layer is not split. |
| `patchBoundaryLayers/cone, walls, tail` | `nLayers 10`, `thicknessRatio 1.2` | Layers on the body. `Allmesh --nLayers` and `--thicknessRatio` change them. |
| `patchBoundaryLayers/fins` | `nLayers 3`, `thicknessRatio 1.2`, `allowDiscontinuity 1` | Layers on the fins. More layers caused tangled cells at the leading edge of the 6 mm fin. `allowDiscontinuity` lets the number of layers change at the fin-body junction. |
| `optimiseLayer` | 1 | Smooths the layers after the extrusion. |
| `untangleLayers` | 1 | Repairs tangled layer cells. |
| `optimisationParameters` | `nSmoothNormals 5`, `maxNumIterations 5`, `featureSizeFactor 0.3`, `reCalculateNormals 1`, `relThicknessTol 0.03` | Parameters of the layer optimization. |

> **CAUTION** cfMesh splits the layer only on the patches that have an
> entry in `patchBoundaryLayers`. With only the global `nLayers 10`, the body
> got no layers, and the faces of the domain got 10. Keep an entry for each
> wall patch.

The layer rule is t₁ = h (r − 1) / (rⁿ − 1), with h the surface cell size, r
the `thicknessRatio` and n the `nLayers`, if the whole layer were one surface
cell thick. `cartesianMesh` makes it thinner: the measured first cell is
0.41 times the rule, with 10 or 12 layers alike. The rule therefore sets the
same t₁ on every mesh, but not its value. Measure it with
`post/firstLayer.py <mesh-dir>`, which gives y (cell centre to wall face, as
`yPlus` uses it) and t₁ = 2y per wall patch, weighted by face area.

| Mesh | `maxCellSize` | h on the body | `nLayers` | `thicknessRatio` | t₁ by the rule | t₁ measured on `walls` |
|---|---|---|---|---|---|---|
| `h0` | 2.2627 | 8.84 mm | 12 | 1.1878 | 0.241 mm | 0.097 mm |
| `h1` | 1.6 | 6.25 mm | 10 | 1.2 | 0.241 mm | 0.100 mm |
| `h2` | 1.1314 | 4.42 mm | 8 | 1.229 | 0.241 mm | |
| `h3` | 0.8 | 3.13 mm | 7 | 1.2 | 0.241 mm | |

So y ≈ 0.05 mm on the body. Where a region refines the wall, the layer is
thinner with it: under `finRoots` (level 11) t₁ is about 0.014 mm, and on
the boattail (`shoulder`, `baseEdge`) about 0.04 mm.

From the flat-plate friction at the aft body (x ≈ 2.8 m), where y+ is
lowest, at the conditions of the sweep:

| M | y+ per mm of y (u_τ/ν) | y+ at the cell centre, y = 0.05 mm |
|---|---|---|
| 0.3 | ~110 | ~5.5 |
| 0.6 | ~230 | ~11.5 |
| 0.8 | ~330 | ~16.5 |
| 1.2 | ~550 | ~27.5 |
| 1.8 | ~950 | ~47.5 |

At M 0.6 to 1.2 the first cell is in the buffer layer, and at M 0.3 at its
lower edge. `nut` is
`nutUSpaldingWallFunction`, which follows Spalding's law through the buffer,
so it holds there. These are estimates: measure y+ on each case.

Procedure to set y+:

1. Make the mesh.
2. Run one case.
3. Read y+ for each patch. Refer to WORKFLOW.md §8.
4. If y+ is not correct, change `nLayers` or `thicknessRatio` and make a new
   mesh.

### 4.5 edgeMeshRefinement

| Entry | Value | Function |
|---|---|---|
| `finEdges/edgeFile` | `"finEdges.obj"` | Leading and trailing edges of the fins. `Allmesh` copies it from the STL directory. |
| `finEdges/additionalRefinementLevels` | 10 | 1.56 mm with `maxCellSize 1.6`. |
| `finEdges/refinementThickness` | 0.010 m | Width of the refined band around each edge. |
| `finRoots/edgeFile` | `"finRoots.obj"` | The root of each fin, LE to TE on the body. `makeRocket.py` writes it, so it follows `--finRoll`. |
| `finRoots/additionalRefinementLevels` | 11 | 0.78 mm with `maxCellSize 1.6`. |
| `finRoots/refinementThickness` | 0.008 m | Width of the band around each root. |

On `h1`, with `finRoots` at level 10 and 0.012 m, `shoulder`, `baseEdge` and `finRoots` added 185 696 cells (1 583 757 → 1 769 453). With `finRoots` at level 11 and 0.008 m, and `surfaceFeatureEdges -angle 60`, `h1` has 2 198 939 cells: non-orthogonality 66.0°, skewness 1.72, 388 faces with interpolation weight below 0.05, 0 negative volumes.

> **CAUTION** `finRoots` refines the body wall as well. The layer thickness follows the surface cell, so at level 11 the first body layer near the root is about 0.014 mm, against 0.10 mm on the rest of the body, and y+ there is about 7 times lower. 98 % of the edges below 20 µm of `h1` are in this band.

The oblique shock starts at the leading edge and the expansion at the
trailing edge.

## 5. Domain

`Allmesh` step 4 sets the domain with six distances in metres from the
bounding box of the rocket:

```
surfaceGenerateBoundingBox surface.stl geometry.stl  17.73 38.415 11.5395 11.5395 11.5395 11.5395
                                                     x-    x+     y-      y+      z-      z+
```

| Direction | Distance from the rocket | Domain limits |
|---|---|---|
| Upstream (−x) | 17.73 m | x = −17.73 m |
| Downstream (+x) | 38.415 m | x = 41.37 m |
| Lateral (±y, ±z) | 11.5395 m | ±11.78 m |

These values are in `Allmesh`, not in a dictionary. They are for the
Aconcagua.

### 5.1 Half mesh

At α = 0 the flow has the symmetry of the fins. By default `Allmesh` meshes
the whole rocket with the fins at 45° (§3.3) and keeps the half z ≥ 0;
`--full` meshes the whole rocket with the fins at 0°.

```
surfaceGenerateBoundingBox surface.stl geometry.stl  17.73 38.415 11.604756 11.604756 11.604756 11.604756
```

The box is symmetric about z = 0, so the octree splits on z = 0 and the
cells of the full mesh do not cross it, except for about 0.1 % near the wall.
`topoSet` selects the cells with the centre at z ≥ 0, and `subsetMesh` keeps
them; the faces of the cut go to `oldInternalFaces`. The smoothing of
`cartesianMesh` leaves the cut up to 1.7 mm off the plane, in cells of 50 to
400 mm, and up to 0.2 mm in the layers. `snapPlane.py` sets z = 0 on every
point of the cut, and `createPatchDictHalf` makes it the patch `symm`, of type
`symmetry`, which `0.orig/` already has. `meshInfo` gets `sector half`,
`symmetryFraction 0.5` and `Aref` = π·rBody²/2, so the coefficients are those
of the full rocket. `Allrun` then permits α but not β.

The half mesh is the full one cut in two. On `h1` the full mesh at 45° has
3 152 544 cells and the half 1 577 819, with skewness 3.58, 156 faces above
70° and 478 faces with interpolation weight below 0.05, half of those of the
full mesh. `cartesianMesh` makes the full mesh, so a half mesh takes the
meshing time and memory of the full one.

> **CAUTION** Do not mesh z ≥ 0 directly with `zMin` on z = 0. `cartesianMesh`
> then treats `zMin` as a domain face that touches the wall: it gives `zMin`
> the layers of the wall, and the wall layers fan out at the corner. On `h1`
> this gave 9 002 faces with interpolation weight below 0.05. The type
> `symmetryPlane` in the `.fms`, `boundaryLayers/excludedPatches` and
> `allowDiscontinuity` did not correct it. At M 1.8 the Mach angle is 33.7°, and the Mach cone from the nose
crosses the lateral boundary at approximately x = 17.6 m. The pressure
boundary condition on `box` must thus let the wave go out. Refer to
SOLVERS.md §9.

## 6. createPatchDict

File: `mesh/system/createPatchDict`. `surfaceGenerateBoundingBox` makes six
patches with fixed names. `createPatch` merges them:

| New patch | Type | Source patches |
|---|---|---|
| `inlet` | `patch` | `xMin` |
| `outlet` | `patch` | `xMax` |
| `box` | `patch` | `yMin`, `yMax`, `zMin`, `zMax` |

`createPatchDictHalf`, for the half mesh, makes the faces of the cut at z = 0
(`oldInternalFaces`, from `subsetMesh`) the patch `symm`, of type `symmetry`.

`pointSync false`. The file does not change the wall patches.

### 6.1 Type of the wall patches

The `.fms` file gives no type to its patches. `cartesianMesh` makes all of
them `wall`. `createPatch` makes `inlet`, `outlet` and `box` as new patches
of type `patch`, so only `cone`, `walls`, `tail` and `fins` stay `wall`.

The wall functions in `0.orig/` and the `yPlus` function object need type
`wall`. After `createPatch`, `Allmesh` reads the type of each patch in
`wallPatches` of `meshInfo` from `constant/polyMesh/boundary`. If a type is
not `wall`, `Allmesh` stops.

## 7. meshInfo

File: `mesh/meshInfo`. `Allmesh` copies it into the mesh directory.
`newCase.sh` copies it to `constant/meshInfo` in the case. It replaces the
`meshInfo` of the template. The case reads only this copy.

| Entry | Value | Written by | Read by | Function |
|---|---|---|---|---|
| `sector` | `full` | user, `Allmesh` (half) | `Allrun` | Part of the rocket in the mesh: `quarter`, `half` or `full`. Sets the permitted flow angles. |
| `symmetryFraction` | 1 | user, `Allmesh` (half) | no script | Fraction of the rocket in the mesh. |
| `fins` | `true` | user | `Allrefine` | `Allrefine` stops if it is not `true`. |
| `nCells` | from `checkMesh` | `Allmesh`, `Allrefine` | no script | Number of cells. |
| `Aref` | 0.0179081 m² | user, `Allmesh` (half) | `controlDictBase` | Reference area: π·rBody² for the full rocket, half of it for a half mesh. |
| `lRef` | 0.151 m | user | `controlDictBase` | Reference length: body diameter. |
| `rBody` | 0.0755 m | user | `Allmesh` → `meshDict` | Body radius. Scales the regions across the axis (§4.3). |
| `lBody` | 2.955 m | user | `Allmesh` → `meshDict` | Rocket length. Scales the regions along the axis (§4.3). |
| `finRootR` | 0.075 m | user | no script | Radius of the fin root. |
| `finTipR` | 0.2355 m | user | `topoSetDict`, `Allrefine` | Radius of the fin tip. |
| `finTipSmear` | 0.004 m | user | no script | Width of the tip chamfer. |
| `finX0` | 2.52953 m | user | `topoSetDict` | x of the fin front. |
| `finX1` | 2.93031 m | user | `topoSetDict` | x of the fin rear. |
| `wallPatches` | `(cone walls tail fins)` | user | `controlDictBase` | Patches in the force integration and in `wallShearStress`. |
| `maxCellSize` | from `meshDict` | `Allmesh` | no script | Record of the mesh. |
| `finLevel` | from `meshDict` | `Allmesh` | no script | Record of the mesh. |
| `nLayers`, `thicknessRatio` | from `meshDict`, patch `walls` | `Allmesh` | no script | Record of the mesh. |
| `meshDir` | path | `Allmesh` | no script | Record of the mesh. |

> **NOTE** `Aref` must agree with `sector`. For a `half` mesh, use half the
> area. For a `quarter` mesh, use a quarter of the area. The
> `case-central/constant/meshInfo` file has the values of a `quarter` mesh.
> `newCase.sh` replaces it, so these values have no effect.

## 8. Other files in mesh/system/

`controlDict`, `fvSchemes` and `fvSolution` are necessary because the
OpenFOAM applications read them at start. They have no effect on the mesh.

## 9. Examine a mesh

1. Read the `Allmesh` summary. There must be no `negative cell volume` and
   no `Failed` check.
2. Examine the patch names and types in `constant/polyMesh/boundary`.
   Refer to §6.1.
3. Open the mesh in ParaView. Look at the cells on the fin edges, on the
   base edge and on the boattail shoulder.
4. Make sure that the number of cells is correct for the case cost.
