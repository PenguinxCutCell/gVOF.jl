"""
    gVOF

A Julia package implementing advanced unsplit geometric Volume of Fluid (VOF)
methods on arbitrary grids, either structured or unstructured with convex or
non-convex cells.

Translated from the Fortran gVOF library by J. Lopez and J. Hernandez (2021).

# References
- J. Lopez and J. Hernandez, gVOF: An open-source package for unsplit geometric
  volume of fluid methods on arbitrary grids, Computer Physics Communications.

# Dependencies
- `ISOAP` – isosurface extraction on arbitrary polyhedra
- `VOFTools` – analytical and geometrical tools for VOF methods
"""
module gVOF

using ISOAP
using ISOAP: niso
using VOFTools
using LinearAlgebra
using Printf

# ── Re-exports from dependencies ──────────────────────────────────────────────
export Polyhedron, Grid, IsoResult, isovtkgrid
export Polyhedron3D

# ── Types ─────────────────────────────────────────────────────────────────────
export VOFGrid, VOFParams

# ── Grid operations ───────────────────────────────────────────────────────────
export vofgrid, compgrid!, neigbcell!, defcell

# ── Tagging & initialisation ─────────────────────────────────────────────────
export taggrid!, initfgrid!, recerr

# ── Reconstruction ────────────────────────────────────────────────────────────
export clcir!, elcir!, llcir!, lsgir!, swir!, lsfir!

# ── Advection ─────────────────────────────────────────────────────────────────
export emfp!, fmfp!, nmfp!, faceflux!, vofadv!

# ── Output ────────────────────────────────────────────────────────────────────
export printgrid, printplic, printvoxel

# ── Test cases ────────────────────────────────────────────────────────────────
export vardef, velgrid!, compdt!
export func3dinic, func3dendc, funcexactvol

# ── Convenience ───────────────────────────────────────────────────────────────
export reconstruct!

# ── Source files ──────────────────────────────────────────────────────────────
include("types.jl")
include("grid.jl")
include("reconstruction.jl")
include("advection.jl")
include("tagging.jl")
include("output.jl")
include("testcases.jl")

end # module gVOF
