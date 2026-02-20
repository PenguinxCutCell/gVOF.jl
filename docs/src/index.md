```@meta
CurrentModule = gVOF
```

# gVOF

## Installation

```julia
using Pkg
Pkg.add("gVOF")
```

## Overview

`gVOF.jl` is a Julia implementation of unsplit geometric VOF methods for
structured and unstructured grids.

Main capabilities:

- Grid preprocessing and geometric metrics
- Interface tagging and volume-fraction initialization
- Multiple PLIC reconstruction methods
- Geometric advection with flux polyhedra
- VTK output helpers for grid/interface/voxel fields

## Quick Example

```julia
using gVOF

params = VOFParams(nx = 32, ny = 32, nz = 32, igrid = 1, icase = 3)
vg = vofgrid(params)
initfgrid!(vg, func3dinic(params.icase); nc = params.nc)
taggrid!(vg; tolfr = params.tolfr)
reconstruct!(vg, params)
```

## Main Sections

- [Types](@ref)
- [Grid](@ref)
- [Tagging and Initialization](@ref)
- [Reconstruction](@ref)
- [Advection](@ref)
- [Output and Test Cases](@ref)
- [Reference](@ref reference)

## Performance

Allocations free and fast execution.

## Reference

Joaquín López, Julio Hernández, gVOF: An open-source package for unsplit geometric volume of fluid methods on arbitrary grids, Computer Physics Communications, Volume 277, 2022, 108400, ISSN 0010-4655, https://doi.org/10.1016/j.cpc.2022.108400.