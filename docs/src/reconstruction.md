```@meta
CurrentModule = gVOF
```

# Reconstruction

ELVIRA and LVIRA are currently implemented as canonical 2D methods for
uniform Cartesian grids with a single z-layer (`nz == 1`). They use the same
3D polyhedron clipping path as the other reconstructions, but the interface
normal is constrained to the x-y plane and the interface position is enforced
exactly in the center cell.

ELVIRA evaluates the six canonical slope candidates derived from the 3×3
stencil column and row sums, enforces the center-cell volume exactly for each
candidate, and selects the one with the smallest stencil mismatch.

LVIRA parameterizes the unit normal by a single angle, enforces the same
center-cell equality constraint, and minimizes the stencil mismatch with a
dependency-free 1D search.

Cells that do not have a full 3×3 stencil fall back locally to a 2D
least-squares gradient reconstruction. Unsupported grids raise a clear error.

```@docs
clcir!
elcir!
elvira!
llcir!
lsgir!
lvira!
swir!
lsfir!
reconstruct!
```
