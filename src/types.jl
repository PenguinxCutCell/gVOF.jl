# ---------------------------------------------------------------------------
#  Types for gVOF.jl
# ---------------------------------------------------------------------------

"""
    VOFParams

Parameters controlling a gVOF simulation, typically read from a `gvofvardef`
configuration file.
"""
@kwdef struct VOFParams
    icase::Int       = 3       # Test case number
    irec::Int        = 1       # Interface reconstruction method (1–6)
    pn::Float64      = 1.5     # Distance weighting parameter
    iw::Int          = 1       # Iso-surface weighting factor (1–6)
    niter::Int       = 4       # Outer iterations for SWIR / LSFIR
    tolir::Float64   = 0.001   # Tolerance for SWIR / LSFIR
    iadv::Int        = 1       # Advection method (1=EMFPA, 2=FMFPA, 3=NMFPA)
    tcaden::Float64  = 0.1     # Print cadence
    iprintplic::Bool = true    # Print PLIC results
    iprintiso::Bool  = false   # Print iso-surface results
    iprintvoxel::Bool = false  # Print voxelisation results
    nc::Int          = 10      # Sub-cells for initialisation / error computation
    tolfr::Float64   = 1e-12   # Volume fraction tolerance
    idt::Int         = 1       # Time-step type (0=fixed, 1=CFL)
    dt::Float64      = 1e-5    # Fixed time step
    cfl::Float64     = 0.5     # CFL number
    igrid::Int       = 1       # Grid type (1=uniform, 2=convex OF, 3=non-convex OF)
    nx::Int          = 32
    ny::Int          = 32
    nz::Int          = 32
    xlast::Float64   = 1.0     # Grid size along X
    ylast::Float64   = 1.0     # Grid size along Y
    zlast::Float64   = 1.0     # Grid size along Z
    timeend::Float64 = 3.0     # Simulation end time
end

"""
    VOFGrid

A grid augmented with all the pre-computed quantities needed by the gVOF
algorithms: face normals/areas, cell centres/volumes, bounding boxes,
neighbour lists, etc.

Built from an `ISOAP.Grid` via [`vofgrid`](@ref).
"""
mutable struct VOFGrid
    # ── base grid from ISOAP ──────────────────────────────────────────────
    grid::Grid

    # ── per-node: cells sharing each node ─────────────────────────────────
    icnode::Vector{Vector{Int}}   # icnode[node] = [cell indices…]

    # ── per-cell: neighbour cells ─────────────────────────────────────────
    ineigb::Vector{Vector{Int}}   # ineigb[cell] = [neighbour indices…]

    # ── per-face quantities ───────────────────────────────────────────────
    aface::Vector{Float64}        # face areas
    cface::Matrix{Float64}        # face centroids   (nface × 3)
    xnface::Vector{Float64}       # face normal x-component
    ynface::Vector{Float64}       # face normal y-component
    znface::Vector{Float64}       # face normal z-component
    nfaceint::Int                 # number of interior faces

    # ── per-cell quantities ───────────────────────────────────────────────
    boxcell::Matrix{Float64}      # bounding boxes  (ncell × 6)
    ccell::Matrix{Float64}        # cell centroids   (ncell × 3)
    vcell::Vector{Float64}        # cell volumes
    sizecmin::Vector{Float64}     # minimum cell sizes (3)

    # ── VOF state ─────────────────────────────────────────────────────────
    fractg::Vector{Float64}       # material volume fraction
    frnod::Vector{Float64}        # nodal volume fraction (interpolated)
    xnormg::Vector{Float64}       # PLIC normal x
    ynormg::Vector{Float64}       # PLIC normal y
    znormg::Vector{Float64}       # PLIC normal z
    rholig::Vector{Float64}       # PLIC plane constant

    # ── tagging arrays ────────────────────────────────────────────────────
    ictag::Vector{Int}            # cell tag (-1, 0, 1)
    iptag::Vector{Int}            # node tag (-1, 0, 1)
    istag::Vector{Int}            # face tag (-1, 0, 1)
    icint::Vector{Int}            # list of interfacial cells
    icadv::Vector{Int}            # list of advection cells
    isflu::Vector{Int}            # list of faces requiring flux computation

    # ── advection arrays ──────────────────────────────────────────────────
    fface::Vector{Float64}        # face liquid flux
    volpol::Vector{Float64}       # flux polyhedron volume
    velface::Matrix{Float64}      # face velocities      (nface × 3)
    velnode::Matrix{Float64}      # node velocities      (npoint × 3)
    velcmax::Vector{Float64}      # max velocity components (3)
    ebound::Vector{Float64}       # boundedness errors (2)
end
