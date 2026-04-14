```@meta
CurrentModule = gVOF
```

# Reconstruction

ELVIRA is currently implemented as a canonical 2D method for
uniform Cartesian grids with a single z-layer (`nz == 1`). They use the same
3D polyhedron clipping path as the other reconstructions, but the interface
normal is constrained to the x-y plane and the interface position is enforced
exactly in the center cell.

ELVIRA evaluates the six canonical slope candidates derived from the 3×3
stencil column and row sums, enforces the center-cell volume exactly for each
candidate, and selects the one with the smallest stencil mismatch.

The 2D LVIRA helper (`lvira!`) parameterizes the unit normal by a single angle,
enforces the same center-cell equality constraint, and minimizes the stencil
mismatch with a dependency-free 1D search.

`lvira3d!` is a full 3D equality-constrained LVIRA formulation for arbitrary
polyhedral grids. It parameterizes the unit normal with two angles, enforces
the center-cell volume exactly, and evaluates the same reconstructed plane
across immediate neighbors. The objective is minimized with a derivative-free
local angular search.

Cells that do not have a full 3×3 stencil fall back locally to a 2D
least-squares gradient reconstruction. Unsupported grids raise a clear error.

```@docs
clcir!
elcir!
elvira!
llcir!
lsgir!
lvira!
lvira3d!
swir!
lsfir!
reconstruct!
```
