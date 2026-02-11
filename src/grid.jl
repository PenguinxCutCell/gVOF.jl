# ---------------------------------------------------------------------------
#  Grid construction and utilities  (translation of DEFGRID, COMPGRID,
#  NEIGBCELL, DEFCELL from gvof.f)
# ---------------------------------------------------------------------------

"""
    vofgrid(params::VOFParams) -> VOFGrid

Build a `VOFGrid` from simulation parameters.  Delegates the raw grid
construction to `ISOAP.constgrid` and then pre-computes every auxiliary
quantity needed by the VOF algorithms.
"""
function vofgrid(params::VOFParams)
    # Build the ISOAP grid
    g = constgrid(params.igrid, params.nx, params.ny, params.nz,
                  params.xlast, params.ylast, params.zlast)
    vg = _allocate_vofgrid(g)
    compgrid!(vg)
    neigbcell!(vg)
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
    icnode = [Int[] for _ in 1:npoint]
    for ic in 1:ncell
        for ip in g.ipcell[ic]
            push!(icnode[ip], ic)
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

    for ic in 1:ncell
        nodes = g.ipcell[ic]
        nip   = length(nodes)

        # Gather local vertex coordinates and build local-index map
        vertp_local = Matrix{Float64}(undef, nip, 3)
        ipg = Vector{Int}(undef, nip)
        xv = yv = zv = 0.0
        xmin = ymin = zmin =  1e16
        xmax = ymax = zmax = -1e16
        for k in 1:nip
            gp = nodes[k]
            ipg[k] = gp
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
        nipv = zeros(Int, nts)
        ipv  = zeros(Int, nts, nip)
        xns  = zeros(nts); yns = zeros(nts); zns = zeros(nts)

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
                local_idx = findfirst(==(gp_face), ipg)
                ipv[is, iv] = local_idx
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
    for ic in 1:g.ncell
        neighbours = Set{Int}()
        for ip in g.ipcell[ic]
            for icn in vg.icnode[ip]
                icn != ic && push!(neighbours, icn)
            end
        end
        vg.ineigb[ic] = collect(neighbours)
    end
    return nothing
end

# -----------------------------------------------------------------------
#  DEFCELL – build a Polyhedron3D for a single grid cell
# -----------------------------------------------------------------------
"""
    defcell(vg::VOFGrid, ic::Int) -> Polyhedron3D

Construct the local `Polyhedron3D` representation of cell `ic`.  Face
winding is reversed for faces where `ic` is the neighbour cell so that
all face normals point outward.
"""
function defcell(vg::VOFGrid, ic::Int)
    g     = vg.grid
    vnode = g.vnode
    nodes = g.ipcell[ic]
    nip   = length(nodes)

    # Build local vertices
    vertp = Matrix{Float64}(undef, nip, 3)
    ipg   = Vector{Int}(undef, nip)
    for k in 1:nip
        gp = nodes[k]
        ipg[k] = gp
        vertp[k,1] = vnode[gp,1]
        vertp[k,2] = vnode[gp,2]
        vertp[k,3] = vnode[gp,3]
    end

    nts  = length(g.iscell[ic])
    nipv = zeros(Int, nts)
    ipv  = zeros(Int, nts, nip)
    xns  = zeros(nts); yns = zeros(nts); zns = zeros(nts)

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
            local_idx = findfirst(==(gp_face), ipg)
            ipv[is, iv] = local_idx
        end
    end

    # Pad to VOFTools layout (use VOFTools defaults for buffer safety)
    nv_max = VOFTools.NV_DEFAULT
    ns_max = VOFTools.NS_DEFAULT
    ipv_pad  = zeros(Int, ns_max, nv_max)
    nipv_pad = zeros(Int, ns_max)
    v_pad    = zeros(nv_max, 3)
    xns_pad  = zeros(ns_max); yns_pad = zeros(ns_max); zns_pad = zeros(ns_max)

    for is in 1:nts
        nipv_pad[is] = nipv[is]
        xns_pad[is]  = xns[is]; yns_pad[is] = yns[is]; zns_pad[is] = zns[is]
        for iv in 1:nipv[is]
            ipv_pad[is, iv] = ipv[is, iv]
        end
    end
    for k in 1:nip
        v_pad[k,1] = vertp[k,1]; v_pad[k,2] = vertp[k,2]; v_pad[k,3] = vertp[k,3]
    end

    return Polyhedron3D(v_pad, ipv_pad, nipv_pad, xns_pad, yns_pad, zns_pad,
                        nts, nip, nip)
end
