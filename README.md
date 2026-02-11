# gVOF.jl

A Julia implementation of the geometrical Volume-of-Fluid (gVOF) method for interface tracking and reconstruction.

**Original Work:** This package is a Julia translation of the Fortran gVOF library by J. Lopez and J. Hernandez (2021). All credits for the original research, algorithms, and Fortran implementation go to the original authors and contributors.

**Reference:**
- J. Lopez and J. Hernandez, "gVOF: An open-source package for unsplit geometric volume of fluid methods on arbitrary grids," *Computer Physics Communications*.

## Features

- Accurate interface reconstruction on arbitrary grids
- Volume fraction computation and advection
- Support for multi-dimensional flows on structured and unstructured grids
- Six reconstruction methods (CLCIR, ELCIR, LLCIR, LSGIR, SWIR, LSFIR)
- Three advection schemes (EMFP, FMFP, NMFP)
- Efficient Julia-based implementation

## Installation

```julia
using Pkg
Pkg.add("gVOF")
```

## Usage

```julia
using gVOF

# Example: Initialize and run gVOF solver
params = VOFParams(nx=32, ny=32, nz=32)
vg = vofgrid(params)
initfgrid!(vg, func3dinic(1))
taggrid!(vg)
reconstruct!(vg, params)
```

