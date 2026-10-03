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
| 2 | `foamDictionary` | `--maxCellSize`, `--finLevel` | `system/meshDict` | Writes the options into the copy. |
| 3 | `sed` | STL | `surface.stl` | Changes the solid names to the patch names. Refer to §3. |
| 4 | `surfaceGenerateBoundingBox` | `surface.stl` | `geometry.stl` | Adds the six faces of the domain. Refer to §5. |
| 5 | `surfaceFeatureEdges -angle 30` | `geometry.stl` | `geometry.fms` | Finds the feature edges. Refer to §2.1. |
| 6 | `cartesianMesh` | `geometry.fms`, `system/meshDict` | `constant/polyMesh` | Makes the mesh. |
| 7 | `createPatch -overwrite` | `system/createPatchDict` | `constant/polyMesh` | Merges the six domain faces into `inlet`, `outlet` and `box`. |
| 8 | `checkMesh -allTopology -allGeometry` | mesh | `log.checkMesh` | Measures the mesh quality. |
| 9 | `foamDictionary` | `log.checkMesh`, `meshDict` | `meshInfo` | Writes `nCells`, `maxCellSize`, `finLevel`, `meshDir`. |
| 10 | `grep` | `log.checkMesh` | screen | Shows the summary. Stops if there are negative cell volumes. |

### 2.1 Function of the .fms file

An STL file contains only triangles and solid names. It does not contain
edges. The `.fms` file contains the triangles, the patches and the feature
edges.

`surfaceFeatureEdges` marks an edge as a feature edge when the angle between
the two adjacent triangles is more than 30°. `cartesianMesh` then moves the
mesh vertices onto these edges.

If an edge is not a feature edge, the mesh cuts it at the cell size. The
edge then becomes round or chamfered. These edges have an effect on the
results:

- The leading and trailing edges and the tips of the fins. They set the
  position of the shocks and the flow separation on the fins.
- The fin-body junction (90°).
- The edge of the base. It sets the separation point and thus the base
  pressure.

> **NOTE** An edge with an angle of less than 30° is not a feature edge. The
> boattail shoulder has the angle of the boattail, which is usually less
> than 30°. A round leading edge made of small facets also has small angles.
> To find the feature edges, run `surfaceFeatureEdges` with a smaller angle
> and look at the result in ParaView.

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

1. Put the STL in `mesh/stl/`.
2. Make sure that the STL agrees with the requirements in §3.
3. Change the distances in `surfaceGenerateBoundingBox` in `Allmesh`. Refer
   to §5.
4. Change the positions of the regions in `objectRefinements` in
   `meshDict`. Refer to §4.3.
5. Change `mesh/meshInfo`. Refer to §7.
6. Run `mesh/Allmesh <mesh-dir> mesh/stl/<file>.stl`.

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
| `fins` | 9 | 3.13 mm | 0.020 m | Fin chord. `Allmesh --finLevel` changes this level. |

### 4.3 objectRefinements

Volume refinement in regions. A `box` has a centre and three lengths. A
`cone` has two points on the axis and a radius at each point.

| Region | Type | Level | Cell size | Extent |
|---|---|---|---|---|
| `nearField` | box | 5 | 50 mm | x from −0.20 to 3.20 m, y and z ±0.50 m |
| `baseWake` | cone | 8 | 6.25 mm | x from 2.94 to 3.30 m, radius 0.08 m |
| `nearWake` | cone | 7 | 12.5 mm | x from 2.90 to 3.80 m, radius 0.14 m |
| `midWake` | cone | 6 | 25 mm | x from 2.85 to 5.00 m, radius 0.26 m |
| `wake` | box | 4 | 100 mm | x from 1.0 to 9.0 m, y and z ±0.80 m |
| `farWake` | box | 2 | 400 mm | x from 1.0 to 27.0 m, y and z ±2.0 m |

The wake regions decrease the cell size one level at a time behind the base.
The positions are for the Aconcagua (body length 2.955 m). For a different
geometry, change them.

### 4.4 boundaryLayers

cfMesh extrudes the layers from the surface mesh. The first cell thickness
thus comes from the surface cell size, `nLayers` and `thicknessRatio`. The
dictionary does not set a y+ target.

| Entry | Value | Function |
|---|---|---|
| `nLayers` | 10 | Number of layers on all patches that have no entry in `patchBoundaryLayers`. |
| `thicknessRatio` | 1.2 | Ratio between the thicknesses of two adjacent layers. |
| `patchBoundaryLayers/cone, walls, tail` | `nLayers 10` | Layers on the body. |
| `patchBoundaryLayers/fins` | `nLayers 3`, `allowDiscontinuity 1` | Layers on the fins. The fin is thin. More layers cause tangled cells at the leading edge. `allowDiscontinuity` lets the number of layers change at the fin-body junction. |
| `optimiseLayer` | 1 | Smooths the layers after the extrusion. |
| `untangleLayers` | 1 | Repairs tangled layer cells. |
| `optimisationParameters` | `nSmoothNormals 5`, `maxNumIterations 5`, `featureSizeFactor 0.3`, `reCalculateNormals 1`, `relThicknessTol 0.03` | Parameters of the layer optimization. |

Procedure to set y+:

1. Make the mesh.
2. Run one case.
3. Read y+ for each patch. Refer to WORKFLOW.md §8.
4. If y+ is not correct, change `nLayers` or `thicknessRatio` and make a new
   mesh.

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
Aconcagua. At M 1.8 the Mach angle is 33.7°, and the Mach cone from the nose
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

`pointSync false`. The file does not change the wall patches.

### 6.1 Type of the wall patches

> **CAUTION** No step in the procedure sets the type `wall` on `cone`,
> `walls`, `tail` and `fins`. They get the type from the `.fms` file. The
> header of `createPatchDict` refers to a script `fixPatchTypes.py`. This
> script is not in the repository.
>
> The wall functions in `0.orig/` (`nutUSpaldingWallFunction`,
> `kLowReWallFunction`, `omegaWallFunction`, `alphatWallFunction`) and the
> `yPlus` function object need patches of type `wall`. On a patch of type
> `patch`, the solver stops with an error.

Procedure to examine the types:

```bash
grep -A2 -E '^ *(cone|walls|tail|fins)$' <mesh-dir>/constant/polyMesh/boundary
```

If the type is `patch`, one possible correction is the cfMesh entry
`renameBoundary` in `meshDict` [not tested in this repository]:

```
renameBoundary
{
    newPatchNames
    {
        cone  { newName cone;  type wall; }
        walls { newName walls; type wall; }
        tail  { newName tail;  type wall; }
        fins  { newName fins;  type wall; }
    }
}
```

## 7. meshInfo

File: `mesh/meshInfo`. `Allmesh` copies it into the mesh directory.
`newCase.sh` copies it to `constant/meshInfo` in the case. It replaces the
`meshInfo` of the template. The case reads only this copy.

| Entry | Value | Written by | Read by | Function |
|---|---|---|---|---|
| `sector` | `full` | user | `Allrun` | Part of the rocket in the mesh: `quarter`, `half` or `full`. Sets the permitted flow angles. |
| `symmetryFraction` | 1 | user | no script | Fraction of the rocket in the mesh. |
| `fins` | `true` | user | `Allrefine` | `Allrefine` stops if it is not `true`. |
| `nCells` | from `checkMesh` | `Allmesh`, `Allrefine` | no script | Number of cells. |
| `Aref` | 0.0179081 m² | user | `controlDictBase` | Reference area: π·rBody² for the full rocket. |
| `lRef` | 0.151 m | user | `controlDictBase` | Reference length: body diameter. |
| `rBody` | 0.0755 m | user | no script | Body radius. |
| `lBody` | 2.955 m | user | no script | Body length. |
| `finRootR` | 0.075 m | user | no script | Radius of the fin root. |
| `finTipR` | 0.2355 m | user | `topoSetDict`, `Allrefine` | Radius of the fin tip. |
| `finTipSmear` | 0.004 m | user | no script | Width of the tip chamfer. |
| `finX0` | 2.52922 m | user | `topoSetDict` | x of the fin front. |
| `finX1` | 2.93 m | user | `topoSetDict` | x of the fin rear. |
| `wallPatches` | `(cone walls tail fins)` | user | `controlDictBase` | Patches in the force integration and in `wallShearStress`. |
| `maxCellSize` | from `meshDict` | `Allmesh` | no script | Record of the mesh. |
| `finLevel` | from `meshDict` | `Allmesh` | no script | Record of the mesh. |
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
