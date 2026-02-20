# ---------------------------------------------------------------------------
#  Grid construction and utilities  (translation of DEFGRID, COMPGRID,
#  NEIGBCELL, DEFCELL from gvof.f)
# ---------------------------------------------------------------------------

@inline function _find_local_index(ipg::Vector{Int}, nip::Int, gp_face::Int)
    @inbounds for k in 1:nip
        ipg[k] == gp_face && return k
    end
    error("Local vertex index not found for global point $gp_face")
end

struct _VOFGridStatic
    grid::Grid
    icnode::Vector{Vector{Int}}
    ineigb::Vector{Vector{Int}}
    aface::Vector{Float64}
    cface::Matrix{Float64}
    xnface::Vector{Float64}
    ynface::Vector{Float64}
    znface::Vector{Float64}
    nfaceint::Int
    boxcell::Matrix{Float64}
    ccell::Matrix{Float64}
    vcell::Vector{Float64}
    sizecmin::Vector{Float64}
end

@inline function _grid_cache_key(params::VOFParams)
    return (params.igrid, params.nx, params.ny, params.nz,
            params.xlast, params.ylast, params.zlast)
end

function _static_vofgrid(params::VOFParams)
    tls = task_local_storage()
    cache_any = get(tls, :_gvof_static_grid_cache, nothing)
    cache = if cache_any === nothing
        c = Dict{Tuple{Int,Int,Int,Int,Float64,Float64,Float64}, _VOFGridStatic}()
        tls[:_gvof_static_grid_cache] = c
        c
    else
        cache_any::Dict{Tuple{Int,Int,Int,Int,Float64,Float64,Float64}, _VOFGridStatic}
    end
    key = _grid_cache_key(params)
    static = get(cache, key, nothing)
    if static === nothing
        g = constgrid(params.igrid, params.nx, params.ny, params.nz,
                      params.xlast, params.ylast, params.zlast)
        vg = _allocate_vofgrid(g)
        compgrid!(vg)
        neigbcell!(vg)
        static = _VOFGridStatic(
            vg.grid,
            vg.icnode,
            vg.ineigb,
            vg.aface,
            vg.cface,
            vg.xnface,
            vg.ynface,
            vg.znface,
            vg.nfaceint,
            vg.boxcell,
            vg.ccell,
            vg.vcell,
            vg.sizecmin,
        )
        cache[key] = static
    end
    return static::_VOFGridStatic
end

"""
    vofgrid(params::VOFParams) -> VOFGrid

Build a `VOFGrid` from simulation parameters.  Delegates the raw grid
construction to `ISOAP.constgrid` and then pre-computes every auxiliary
quantity needed by the VOF algorithms.
"""
function vofgrid(params::VOFParams)
    static = _static_vofgrid(params)
    return VOFGrid(
        # base grid
        static.grid,
        # per-node
        static.icnode,
        # per-cell neighbours
        static.ineigb,
        # per-face
        static.aface,
        static.cface,
        static.xnface,
        static.ynface,
        static.znface,
        static.nfaceint,
        # per-cell
        static.boxcell,
        static.ccell,
        static.vcell,
        static.sizecmin,
        # VOF state
        zeros(static.grid.ncell),
        zeros(static.grid.npoint),
        zeros(static.grid.ncell),
        zeros(static.grid.ncell),
        zeros(static.grid.ncell),
        zeros(static.grid.ncell),
        # tagging
        zeros(Int, static.grid.ncell),
        zeros(Int, static.grid.npoint),
        zeros(Int, static.grid.nface),
        Int[],
        Int[],
        Int[],
        # advection
        zeros(static.grid.nface),
        zeros(static.grid.nface),
        zeros(static.grid.nface, 3),
        zeros(static.grid.npoint, 3),
        fill(1e-16, 3),
        zeros(2),
    )
end

"""
    vofgrid!(vg::VOFGrid, params::VOFParams) -> VOFGrid

Reset an existing `VOFGrid` in-place for `params`, reusing all mutable
runtime arrays. This is allocation-free when `vg` already matches the
requested grid dimensions.
"""
function vofgrid!(vg::VOFGrid, params::VOFParams)
    static = _static_vofgrid(params)
    g = static.grid

    if vg.grid.ncell != g.ncell || vg.grid.nface != g.nface || vg.grid.npoint != g.npoint
        error("vofgrid! requires a preallocated VOFGrid with matching dimensions")
    end

    # Static geometry/topology
    vg.grid = static.grid
    vg.icnode = static.icnode
    vg.ineigb = static.ineigb
    vg.aface = static.aface
    vg.cface = static.cface
    vg.xnface = static.xnface
    vg.ynface = static.ynface
    vg.znface = static.znface
    vg.nfaceint = static.nfaceint
    vg.boxcell = static.boxcell
    vg.ccell = static.ccell
    vg.vcell = static.vcell
    vg.sizecmin = static.sizecmin

    # Runtime state
    fill!(vg.fractg, 0.0)
    fill!(vg.frnod, 0.0)
    fill!(vg.xnormg, 0.0)
    fill!(vg.ynormg, 0.0)
    fill!(vg.znormg, 0.0)
    fill!(vg.rholig, 0.0)

    fill!(vg.ictag, 0)
    fill!(vg.iptag, 0)
    fill!(vg.istag, 0)
    empty!(vg.icint)
    empty!(vg.icadv)
    empty!(vg.isflu)

    fill!(vg.fface, 0.0)
    fill!(vg.volpol, 0.0)
    fill!(vg.velface, 0.0)
    fill!(vg.velnode, 0.0)
    vg.velcmax[1] = 1e-16
    vg.velcmax[2] = 1e-16
    vg.velcmax[3] = 1e-16
    fill!(vg.ebound, 0.0)

    return vg
end

"""
    _allocate_vofgrid(g::Grid) -> VOFGrid

Allocate a `VOFGrid` wrapping `g` with zeroed arrays.
"""
function _allocate_vofgrid(g::Grid)
    ncell  = g.ncell
    nface  = g.nface
    npoint = g.npoint

    # Per-node: cells sharing each node  (built from ipcell)
    counts = zeros(Int, npoint)
    for ic in 1:ncell
        for ip in g.ipcell[ic]
            counts[ip] += 1
        end
    end
    icnode = Vector{Vector{Int}}(undef, npoint)
    for ip in 1:npoint
        icnode[ip] = Vector{Int}(undef, counts[ip])
        counts[ip] = 0
    end
    for ic in 1:ncell
        for ip in g.ipcell[ic]
            counts[ip] += 1
            icnode[ip][counts[ip]] = ic
        end
    end

    VOFGrid(
        # base grid
        g,
        # per-node
        icnode,
        # per-cell neighbours (empty, filled by neigbcell!)
        [Int[] for _ in 1:ncell],
        # per-face
        zeros(nface),                         # aface
        zeros(nface, 3),                      # cface
        zeros(nface),                         # xnface
        zeros(nface),                         # ynface
        zeros(nface),                         # znface
        0,                                    # nfaceint
        # per-cell
        zeros(ncell, 6),                      # boxcell
        zeros(ncell, 3),                      # ccell
        zeros(ncell),                         # vcell
        fill(1e16, 3),                        # sizecmin
        # VOF state
        zeros(ncell),                         # fractg
        zeros(npoint),                        # frnod
        zeros(ncell),                         # xnormg
        zeros(ncell),                         # ynormg
        zeros(ncell),                         # znormg
        zeros(ncell),                         # rholig
        # tagging
        zeros(Int, ncell),                    # ictag
        zeros(Int, npoint),                   # iptag
        zeros(Int, nface),                    # istag
        Int[],                                # icint
        Int[],                                # icadv
        Int[],                                # isflu
        # advection
        zeros(nface),                         # fface
        zeros(nface),                         # volpol
        zeros(nface, 3),                      # velface
        zeros(npoint, 3),                     # velnode
        fill(1e-16, 3),                       # velcmax
        zeros(2),                             # ebound
    )
end

# -----------------------------------------------------------------------
#  COMPGRID – compute face normals / areas, cell volumes / centroids /
#             bounding boxes
# -----------------------------------------------------------------------
"""
    compgrid!(vg::VOFGrid)

Compute face normals, face areas, face & cell centroids, cell volumes,
bounding boxes and minimum cell sizes.  Mirrors Fortran `COMPGRID`.
"""
function compgrid!(vg::VOFGrid)
    g = vg.grid
    ncell  = g.ncell
    nface  = g.nface
    vnode  = g.vnode

    # ── Face normal & area computation ────────────────────────────────
    for iface in 1:nface
        ipf = g.ipface[iface]
        nip = length(ipf)
        ip1, ip2, ip3 = ipf[1], ipf[2], ipf[3]

        v1 = (vnode[ip2,1]-vnode[ip1,1], vnode[ip2,2]-vnode[ip1,2], vnode[ip2,3]-vnode[ip1,3])
        v2 = (vnode[ip3,1]-vnode[ip2,1], vnode[ip3,2]-vnode[ip2,2], vnode[ip3,3]-vnode[ip2,3])
        xn = v1[2]*v2[3] - v1[3]*v2[2]
        yn = v1[3]*v2[1] - v1[1]*v2[3]
        zn = v1[1]*v2[2] - v1[2]*v2[1]
        d  = sqrt(xn^2 + yn^2 + zn^2)
        vg.xnface[iface] = xn / d
        vg.ynface[iface] = yn / d
        vg.znface[iface] = zn / d

        # Face centroid and area (cross-product formula)
        sx = sy = sz = 0.0
        cx = cy = cz = 0.0
        for i in 1:nip
            ip  = ipf[i]
            jp  = ipf[mod1(i+1, nip)]
            cx += vnode[ip,1]; cy += vnode[ip,2]; cz += vnode[ip,3]
            sx += vnode[ip,2]*vnode[jp,3] - vnode[ip,3]*vnode[jp,2]
            sy += vnode[ip,3]*vnode[jp,1] - vnode[ip,1]*vnode[jp,3]
            sz += vnode[ip,1]*vnode[jp,2] - vnode[ip,2]*vnode[jp,1]
        end
        vg.cface[iface,1] = cx / nip
        vg.cface[iface,2] = cy / nip
        vg.cface[iface,3] = cz / nip
        vg.aface[iface] = (sx*vg.xnface[iface] + sy*vg.ynface[iface] +
                            sz*vg.znface[iface]) / 2.0
    end

    # Count interior faces
    nfint = 0
    for iface in 1:nface
        if g.icface[iface,2] != 0
            nfint += 1
        end
    end
    vg.nfaceint = nfint

    # ── Cell volume, centroid, bounding-box ───────────────────────────
    vt = 0.0
    sizecmin = fill(1e16, 3)
    max_nip = maximum(length, g.ipcell)
    max_nts = maximum(length, g.iscell)
    vertp_local = Matrix{Float64}(undef, max_nip, 3)
    ipg = Vector{Int}(undef, max_nip)
    local_idx = zeros(Int, g.npoint)
    nipv = Vector{Int}(undef, max_nts)
    ipv = Matrix{Int}(undef, max_nts, max_nip)
    xns = Vector{Float64}(undef, max_nts)
    yns = Vector{Float64}(undef, max_nts)
    zns = Vector{Float64}(undef, max_nts)

    for ic in 1:ncell
        nodes = g.ipcell[ic]
        nip   = length(nodes)

        # Gather local vertex coordinates and build local-index map
        xv = yv = zv = 0.0
        xmin = ymin = zmin =  1e16
        xmax = ymax = zmax = -1e16
        for k in 1:nip
            gp = nodes[k]
            ipg[k] = gp
            local_idx[gp] = k
            x, y, z = vnode[gp,1], vnode[gp,2], vnode[gp,3]
            vertp_local[k,1] = x; vertp_local[k,2] = y; vertp_local[k,3] = z
            xv += x; yv += y; zv += z
            xmin = min(xmin, x); xmax = max(xmax, x)
            ymin = min(ymin, y); ymax = max(ymax, y)
            zmin = min(zmin, z); zmax = max(zmax, z)
        end
        vg.boxcell[ic,1] = xmin; vg.boxcell[ic,2] = xmax
        vg.boxcell[ic,3] = ymin; vg.boxcell[ic,4] = ymax
        vg.boxcell[ic,5] = zmin; vg.boxcell[ic,6] = zmax
        sizecmin[1] = min(sizecmin[1], xmax - xmin)
        sizecmin[2] = min(sizecmin[2], ymax - ymin)
        sizecmin[3] = min(sizecmin[3], zmax - zmin)
        vg.ccell[ic,1] = xv / nip
        vg.ccell[ic,2] = yv / nip
        vg.ccell[ic,3] = zv / nip

        # Build local Polyhedron3D and compute volume via toolv3d
        nts  = length(g.iscell[ic])
        for is in 1:nts
            iface = g.iscell[ic][is]
            ipf   = g.ipface[iface]
            nipv[is] = length(ipf)

            if ic == g.icface[iface,1]
                isign = 1; order = 1:nipv[is]
            else
                isign = -1; order = nipv[is]:-1:1
            end
            xns[is] = isign * vg.xnface[iface]
            yns[is] = isign * vg.ynface[iface]
            zns[is] = isign * vg.znface[iface]

            for (iv, idx) in enumerate(order)
                gp_face = ipf[idx]
                ipv[is, iv] = local_idx[gp_face]
            end
        end

        vol = toolv3d(ipv, nipv, nts, vertp_local, xns, yns, zns)
        vg.vcell[ic] = vol
        vt += vol
    end

    vg.sizecmin .= sizecmin
    @info "Total grid volume: $vt"
    return nothing
end

# -----------------------------------------------------------------------
#  NEIGBCELL – find neighbours (cells sharing at least one node)
# -----------------------------------------------------------------------
"""
    neigbcell!(vg::VOFGrid)

Populate `vg.ineigb` with neighbour-cell lists (cells sharing ≥ 1 node).
"""
function neigbcell!(vg::VOFGrid)
    g = vg.grid
    seen = falses(g.ncell)
    buf = Vector{Int}(undef, g.ncell)
    for ic in 1:g.ncell
        nnei = 0
        for ip in g.ipcell[ic]
            for icn in vg.icnode[ip]
                if icn != ic && !seen[icn]
                    nnei += 1
                    buf[nnei] = icn
                    seen[icn] = true
                end
            end
        end
        neigh = vg.ineigb[ic]
        resize!(neigh, nnei)
        for i in 1:nnei
            icn = buf[i]
            neigh[i] = icn
            seen[icn] = false
        end
    end
    return nothing
end

# -----------------------------------------------------------------------
#  DEFCELL – build a Polyhedron3D for a single grid cell
# -----------------------------------------------------------------------
mutable struct _DefcellWork
    poly::Polyhedron3D
    ipg::Vector{Int}
end

function _DefcellWork(vg::VOFGrid)
    ns_max = VOFTools.NS_DEFAULT
    nv_max = VOFTools.NV_DEFAULT
    max_nip = maximum(length, vg.grid.ipcell)
    return _DefcellWork(
        Polyhedron3D(
            zeros(nv_max, 3),
            zeros(Int, ns_max, nv_max),
            zeros(Int, ns_max),
            zeros(ns_max),
            zeros(ns_max),
            zeros(ns_max),
            0, 0, 0,
        ),
        zeros(Int, max_nip),
    )
end

function _get_defcell_work(vg::VOFGrid)
    tls = task_local_storage()
    cache_any = get(tls, :_gvof_defcell_work, nothing)
    cache = if cache_any === nothing
        c = Dict{UInt, _DefcellWork}()
        tls[:_gvof_defcell_work] = c
        c
    else
        cache_any::Dict{UInt, _DefcellWork}
    end
    key = objectid(vg)
    work = get(cache, key, nothing)
    if work === nothing
        work = _DefcellWork(vg)
        cache[key] = work
    end
    return work::_DefcellWork
end

function defcell!(poly::Polyhedron3D, ipg::Vector{Int}, vg::VOFGrid, ic::Int)
    g     = vg.grid
    vnode = g.vnode
    nodes = g.ipcell[ic]
    nip   = length(nodes)

    # Build local vertices
    for k in 1:nip
        gp = nodes[k]
        ipg[k] = gp
        poly.vertp[k,1] = vnode[gp,1]
        poly.vertp[k,2] = vnode[gp,2]
        poly.vertp[k,3] = vnode[gp,3]
    end

    nts  = length(g.iscell[ic])
    poly.nts = nts
    poly.ntp = nip
    poly.ntv = nip

    for is in 1:nts
        iface = g.iscell[ic][is]
        ipf   = g.ipface[iface]
        nfacev = length(ipf)
        poly.nipv[is] = nfacev

        if ic == g.icface[iface,1]
            isign = 1
            poly.xns[is] = isign * vg.xnface[iface]
            poly.yns[is] = isign * vg.ynface[iface]
            poly.zns[is] = isign * vg.znface[iface]
            for iv in 1:nfacev
                gp_face = ipf[iv]
                poly.ipv[is, iv] = _find_local_index(ipg, nip, gp_face)
            end
        else
            isign = -1
            poly.xns[is] = isign * vg.xnface[iface]
            poly.yns[is] = isign * vg.ynface[iface]
            poly.zns[is] = isign * vg.znface[iface]
            for iv in 1:nfacev
                gp_face = ipf[nfacev - iv + 1]
                poly.ipv[is, iv] = _find_local_index(ipg, nip, gp_face)
            end
        end
    end
    return poly
end

"""
    defcell(vg::VOFGrid, ic::Int) -> Polyhedron3D

Construct the local `Polyhedron3D` representation of cell `ic`.  Face
winding is reversed for faces where `ic` is the neighbour cell so that
all face normals point outward.
"""
function defcell(vg::VOFGrid, ic::Int)
    work = _get_defcell_work(vg)
    defcell!(work.poly, work.ipg, vg, ic)
    return cppol3d(work.poly)
end
