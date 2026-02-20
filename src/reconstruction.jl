# ---------------------------------------------------------------------------
#  PLIC interface reconstruction algorithms
#  (CLCIR, ELCIR, LLCIR, LSGIR, SWIR, LSFIR, ISREC, LSREC, SWREC, LSFREC)
# ---------------------------------------------------------------------------

# ═══════════════════════════════════════════════════════════════════════════
#  Helper: ISREC – interface orientation from extracted iso-triangles
# ═══════════════════════════════════════════════════════════════════════════
"""
    isrec(iw, vertc) -> (xn, yn, zn)

Compute the unit-length normal to an iso-surface polygon with centre at
`vertc[end]` and edge mid-points at `vertc[1:end-1]`.

`iw` selects the weighting scheme (1–6, see gVOF documentation).
"""
function isrec(iw::Int, vertc::Vector{NTuple{3,Float64}})
    n   = length(vertc) - 1
    ctr = vertc[end]
    sxn = syn = szn = 0.0

    for i in 1:n
        i2 = (i == n) ? 1 : i + 1
        xv1 = vertc[i][1] - ctr[1]; yv1 = vertc[i][2] - ctr[2]; zv1 = vertc[i][3] - ctr[3]
        xv2 = vertc[i2][1] - ctr[1]; yv2 = vertc[i2][2] - ctr[2]; zv2 = vertc[i2][3] - ctr[3]

        xm = yv1*zv2 - zv1*yv2
        ym = zv1*xv2 - xv1*zv2
        zm = xv1*yv2 - yv1*xv2

        amod1 = sqrt(xv1^2 + yv1^2 + zv1^2)
        amod2 = sqrt(xv2^2 + yv2^2 + zv2^2)
        (amod1 ≤ 1e-16 || amod2 ≤ 1e-16) && continue
        amod  = sqrt(xm^2 + ym^2 + zm^2)
        amod == 0 && continue

        xv1 /= amod1; yv1 /= amod1; zv1 /= amod1
        xv2 /= amod2; yv2 /= amod2; zv2 /= amod2

        w = if iw == 1        # Max weight
                sinα = amod / (amod1 * amod2)
                sinα / (amod1 * amod2)
            elseif iw == 2     # Max + triangle size
                sinα = amod / (amod1 * amod2)
                sinα * amod / (amod1 * amod2)
            elseif iw == 3     # Triangle area
                0.5 * amod
            elseif iw == 4     # Triangle angle
                acos(clamp(xv1*xv2 + yv1*yv2 + zv1*zv2, -1.0, 1.0))
            elseif iw == 5     # Modified triangle angle
                angle = acos(clamp(xv1*xv2 + yv1*yv2 + zv1*zv2, -1.0, 1.0))
                angle < π/2 ? angle : π - angle
            else               # Unweighted
                1.0
            end

        xm /= amod; ym /= amod; zm /= amod
        sxn += w * xm; syn += w * ym; szn += w * zm
    end

    d = sqrt(sxn^2 + syn^2 + szn^2)
    if d != 0
        return (-sxn/d, -syn/d, -szn/d)
    else
        return (1.0, 0.0, 0.0)
    end
end

# ═══════════════════════════════════════════════════════════════════════════
#  Helper: LSREC – weighted least-squares gradient
# ═══════════════════════════════════════════════════════════════════════════
"""
    lsrec(f, vertcx, vertcy, vertcz, pn) -> (xn, yn, zn)

Weighted least-squares estimate of the interface normal.  The last element
of each array is the reference cell.
"""
function lsrec(f::AbstractVector{Float64},
               vertcx::AbstractVector{Float64},
               vertcy::AbstractVector{Float64},
               vertcz::AbstractVector{Float64},
               pn::Float64)
    n = length(f) - 1
    xr, yr, zr, fr = vertcx[n+1], vertcy[n+1], vertcz[n+1], f[n+1]
    ata11 = ata12 = ata13 = 0.0
    ata22 = ata23 = 0.0
    ata33 = 0.0
    atb1 = atb2 = atb3 = 0.0

    for i in 1:n
        d = sqrt((vertcx[i]-xr)^2 + (vertcy[i]-yr)^2 + (vertcz[i]-zr)^2)
        w = 1.0 / d^pn
        a1 = w * (vertcx[i] - xr)
        a2 = w * (vertcy[i] - yr)
        a3 = w * (vertcz[i] - zr)
        b  = w * (f[i] - fr)

        ata11 += a1 * a1
        ata12 += a1 * a2
        ata13 += a1 * a3
        ata22 += a2 * a2
        ata23 += a2 * a3
        ata33 += a3 * a3

        atb1 += a1 * b
        atb2 += a2 * b
        atb3 += a3 * b
    end

    det = ata11*(ata22*ata33 - ata23*ata23) -
          ata12*(ata12*ata33 - ata13*ata23) +
          ata13*(ata12*ata23 - ata13*ata22)
    abs(det) ≤ 1e-30 && return (1.0, 0.0, 0.0)

    inv11 =  (ata22*ata33 - ata23*ata23) / det
    inv12 = -(ata12*ata33 - ata13*ata23) / det
    inv13 =  (ata12*ata23 - ata13*ata22) / det
    inv22 =  (ata11*ata33 - ata13*ata13) / det
    inv23 = -(ata11*ata23 - ata12*ata13) / det
    inv33 =  (ata11*ata22 - ata12*ata12) / det

    xn = inv11*atb1 + inv12*atb2 + inv13*atb3
    yn = inv12*atb1 + inv22*atb2 + inv23*atb3
    zn = inv13*atb1 + inv23*atb2 + inv33*atb3

    if !(isfinite(xn) && isfinite(yn) && isfinite(zn))
        return (1.0, 0.0, 0.0)
    end

    d = sqrt(xn^2 + yn^2 + zn^2)
    if d > 0
        return (xn/d, yn/d, zn/d)
    else
        return (1.0, 0.0, 0.0)
    end
end

# ═══════════════════════════════════════════════════════════════════════════
#  Helper: LSFREC – least-squares plane fitting from PLIC centroids
# ═══════════════════════════════════════════════════════════════════════════
"""
    lsfrec(vertcx, vertcy, vertcz, vnx, vny, vnz) -> (xn, yn, zn)

Least-squares fit reconstruction.  Fits a plane through the PLIC centroid
positions using their normals as constraints.  The last element of each
array is the reference cell.
"""
function lsfrec(vertcx::AbstractVector{Float64},
                vertcy::AbstractVector{Float64},
                vertcz::AbstractVector{Float64},
                vnx::AbstractVector{Float64},
                vny::AbstractVector{Float64},
                vnz::AbstractVector{Float64})
    nplus1 = length(vertcx)
    n = nplus1 - 1
    xr, yr, zr = vertcx[n+1], vertcy[n+1], vertcz[n+1]
    xn0, yn0, zn0 = vnx[n+1], vny[n+1], vnz[n+1]

    # Accumulate weighted moment sums
    sx2 = sxy = sxz = sy2 = syz = sz2 = 0.0
    ncells = 0
    for i in 1:n
        α = acos(clamp(xn0*vnx[i] + yn0*vny[i] + zn0*vnz[i], -1.0, 1.0))
        α < 0.785398 || continue  # < 45°
        xv = vertcx[i] - xr; yv = vertcy[i] - yr; zv = vertcz[i] - zr
        d = max(1e-20, sqrt(xv^2 + yv^2 + zv^2))
        w = 1.0 / d^2.5
        sx2 += xv*xv*w; sxy += xv*yv*w; sxz += xv*zv*w
        sy2 += yv*yv*w; syz += yv*zv*w; sz2 += zv*zv*w
        ncells += 1
    end
    ncells == 0 && return (xn0, yn0, zn0)

    xn, yn, zn = xn0, yn0, zn0
    # Pick the dominant component and solve the 2×2 system
    ax, ay, az = abs(xn0), abs(yn0), abs(zn0)
    if ax >= ay && ax >= az          # X dominant
        det = sy2*sz2 - syz^2
        det == 0 && return (xn0, yn0, zn0)
        yn = (syz*sxz - sz2*sxy) / det
        zn = (syz*sxy - sy2*sxz) / det
        xn = xn0 < 0 ? (yn = -yn; zn = -zn; -1.0) : 1.0
    elseif ay >= ax && ay >= az      # Y dominant
        det = sx2*sz2 - sxz^2
        det == 0 && return (xn0, yn0, zn0)
        xn = (sxz*syz - sz2*sxy) / det
        zn = (sxz*sxy - sx2*syz) / det
        yn = yn0 < 0 ? (xn = -xn; zn = -zn; -1.0) : 1.0
    else                              # Z dominant
        det = sx2*sy2 - sxy^2
        det == 0 && return (xn0, yn0, zn0)
        xn = (sxy*syz - sy2*sxz) / det
        yn = (sxy*sxz - sx2*syz) / det
        zn = zn0 < 0 ? (xn = -xn; yn = -yn; -1.0) : 1.0
    end

    d = sqrt(xn^2 + yn^2 + zn^2)
    d > 0 ? (xn/d, yn/d, zn/d) : (xn0, yn0, zn0)
end

# ═══════════════════════════════════════════════════════════════════════════
#  Helper: volume enforcement dispatch
# ═══════════════════════════════════════════════════════════════════════════
"""
    _enforce_volume(vg, ic, poly, xnc, ync, znc) -> c

Choose `enforv3dsz` for uniform grids and `enforv3d` for general grids.
"""
function _enforce_volume(vg::VOFGrid, ic::Int, poly::Polyhedron3D,
                         xnc::Float64, ync::Float64, znc::Float64;
                         igrid::Int=1)
    v  = vg.fractg[ic] * vg.vcell[ic]
    vt = vg.vcell[ic]
    if igrid == 1
        dx = vg.boxcell[ic,2] - vg.boxcell[ic,1]
        dy = vg.boxcell[ic,4] - vg.boxcell[ic,3]
        dz = vg.boxcell[ic,6] - vg.boxcell[ic,5]
        return enforv3dsz(dx, dy, dz, v, poly.vertp, xnc, ync, znc)
    else
        return enforv3d(poly, v, vt, xnc, ync, znc)
    end
end

# ═══════════════════════════════════════════════════════════════════════════
#  Helper: compute PLIC centroid after inte3d!
# ═══════════════════════════════════════════════════════════════════════════
function _plic_centroid(poly::Polyhedron3D, nts_before::Int)
    if poly.nts == nts_before + 1
        is = poly.nts
        nv = poly.nipv[is]
        cx = cy = cz = 0.0
        for iv in 1:nv
            j = poly.ipv[is, iv]
            cx += poly.vertp[j,1]; cy += poly.vertp[j,2]; cz += poly.vertp[j,3]
        end
        return (cx/nv, cy/nv, cz/nv)
    end
    return nothing
end

mutable struct _ReconWork
    poly::Polyhedron3D
    poly1::Polyhedron3D
    poly2::Polyhedron3D
    ipg::Vector{Int}
    fv::Vector{Float64}
    vx::Vector{Float64}
    vy::Vector{Float64}
    vz::Vector{Float64}
    vcx::Vector{Float64}
    vcy::Vector{Float64}
    vcz::Vector{Float64}
    nnx::Vector{Float64}
    nny::Vector{Float64}
    nnz::Vector{Float64}
    vertc::Vector{NTuple{3,Float64}}
    phi_v::Vector{Float64}
    iplc::Vector{Int}
    nlc_v::Vector{Int}
    ineigblc::Vector{Vector{Int}}
    xlc::Vector{Vector{Float64}}
    ylc::Vector{Vector{Float64}}
    zlc::Vector{Vector{Float64}}
    cplic::Matrix{Float64}
    cplic2::Matrix{Float64}
    mark::Vector{Int}
    mark2::Vector{Int}
    xnorm2::Vector{Float64}
    ynorm2::Vector{Float64}
    znorm2::Vector{Float64}
    rholi2::Vector{Float64}
    alpha::Vector{Float64}
    iso_cells::Vector{Polyhedron}
    iso_ws::IsoapWorkspace
end

function _ReconWork(vg::VOFGrid)
    ns_max = VOFTools.NS_DEFAULT
    nv_max = VOFTools.NV_DEFAULT
    max_neigs = maximum(length, vg.ineigb)
    max_nip = maximum(length, vg.grid.ipcell)
    ncell = vg.grid.ncell
    max_nts = maximum(length, vg.grid.iscell)
    max_facev = maximum(length, vg.grid.ipface)
    poly_builder() = Polyhedron3D(
        zeros(nv_max, 3),
        zeros(Int, ns_max, nv_max),
        zeros(Int, ns_max),
        zeros(ns_max),
        zeros(ns_max),
        zeros(ns_max),
        0, 0, 0,
    )
    iso_cells = [cellgrid(vg.grid, ic) for ic in 1:ncell]
    dummy_ipv = [collect(1:max_facev) for _ in 1:max_nts]
    dummy_vertp = fill((0.0, 0.0, 0.0), max_nip)
    iso_ws = IsoapWorkspace(Polyhedron(dummy_ipv, dummy_vertp))
    return _ReconWork(
        poly_builder(),
        poly_builder(),
        poly_builder(),
        zeros(Int, max_nip),
        zeros(max_neigs + 1),
        zeros(max_neigs + 1),
        zeros(max_neigs + 1),
        zeros(max_neigs + 1),
        zeros(max_neigs + 1),
        zeros(max_neigs + 1),
        zeros(max_neigs + 1),
        zeros(max_neigs + 1),
        zeros(max_neigs + 1),
        zeros(max_neigs + 1),
        NTuple{3,Float64}[],
        zeros(max_nip),
        zeros(Int, ncell),
        zeros(Int, ncell),
        [Int[] for _ in 1:ncell],
        [Float64[] for _ in 1:ncell],
        [Float64[] for _ in 1:ncell],
        [Float64[] for _ in 1:ncell],
        zeros(ncell, 3),
        zeros(ncell, 3),
        zeros(Int, ncell),
        zeros(Int, ncell),
        zeros(ncell),
        zeros(ncell),
        zeros(ncell),
        zeros(ncell),
        zeros(ncell),
        iso_cells,
        iso_ws,
    )
end

@inline function _ensure_lc_storage!(arr::Vector{Float64}, n::Int)
    resize!(arr, n)
    return arr
end

@inline function _ensure_lc_storage!(arr::Vector{Int}, n::Int)
    resize!(arr, n)
    return arr
end

@inline function _load_vertc!(vertc::Vector{NTuple{3,Float64}},
                              xarr::Vector{Float64},
                              yarr::Vector{Float64},
                              zarr::Vector{Float64},
                              n::Int)
    resize!(vertc, n)
    @inbounds for i in 1:n
        vertc[i] = (xarr[i], yarr[i], zarr[i])
    end
    return vertc
end

function _get_recon_work(vg::VOFGrid)
    tls = task_local_storage()
    cache_any = get(tls, :_gvof_recon_work, nothing)
    cache = if cache_any === nothing
        c = Dict{UInt, _ReconWork}()
        tls[:_gvof_recon_work] = c
        c
    else
        cache_any::Dict{UInt, _ReconWork}
    end
    key = objectid(vg)
    work = get(cache, key, nothing)
    if work === nothing
        work = _ReconWork(vg)
        cache[key] = work
    end
    return work::_ReconWork
end

# ═══════════════════════════════════════════════════════════════════════════
#  1.  LLCIR – Local Level-Contour Interface Reconstruction (irec=3)
# ═══════════════════════════════════════════════════════════════════════════
"""
    llcir!(vg::VOFGrid; iw=1, pn=1.5, igrid=1)

PLIC reconstruction using the local level-contour iso-surface (LLCIR).
"""
function llcir!(vg::VOFGrid; iw::Int=1, pn::Float64=1.5, igrid::Int=1)
    g     = vg.grid
    icint = vg.icint
    work  = _get_recon_work(vg)
    iplc  = work.iplc
    nlc_v = work.nlc_v
    xlc   = work.xlc
    ylc   = work.ylc
    zlc   = work.zlc
    iso_ws = work.iso_ws
    iso_cells = work.iso_cells
    fill!(iplc, 0)
    fill!(nlc_v, 0)

    # Local iso-surface construction
    phivlc = 0.5
    for ind in eachindex(icint)
        ic = icint[ind]
        fnmax = 0.0; fnmin = 1.0
        nodes = g.ipcell[ic]
        nnodes = length(nodes)
        phi_v = @view(work.phi_v[1:nnodes])
        @inbounds for i in 1:nnodes
            phi_v[i] = vg.frnod[nodes[i]]
        end
        for p in phi_v
            fnmax = max(fnmax, p); fnmin = min(fnmin, p)
        end
        fnmax ≤ 0.5 && continue
        fnmin ≥ 0.5 && continue

        isoap!(iso_ws, iso_cells[ic], phi_v, phivlc)
        iso_ws.npoly == 1 || continue

        iplc[ic] = ic
        n_lc = iso_ws.poly_lens[1]
        nlc_v[ic] = n_lc
        xarr = _ensure_lc_storage!(xlc[ic], n_lc + 1)
        yarr = _ensure_lc_storage!(ylc[ic], n_lc + 1)
        zarr = _ensure_lc_storage!(zlc[ic], n_lc + 1)
        xsum = ysum = zsum = 0.0
        for i in 1:n_lc
            ip = iso_ws.polys[1, i]
            v = iso_ws.vertiso[ip]
            xarr[i] = v[1]; yarr[i] = v[2]; zarr[i] = v[3]
            xsum += v[1]; ysum += v[2]; zsum += v[3]
        end
        xarr[end] = xsum / n_lc; yarr[end] = ysum / n_lc; zarr[end] = zsum / n_lc
    end

    # Volume enforcement
    for ind in eachindex(icint)
        ic = icint[ind]
        poly = defcell!(work.poly, work.ipg, vg, ic)

        if iplc[ic] != 0
            n_lc = nlc_v[ic]
            vertc = _load_vertc!(work.vertc, xlc[ic], ylc[ic], zlc[ic], n_lc + 1)
            xn, yn, zn = isrec(iw, vertc)
        else
            neigs = vg.ineigb[ic]
            nn    = length(neigs)
            fv = work.fv
            vx = work.vx
            vy = work.vy
            vz = work.vz
            for ii in 1:nn
                ic2 = neigs[ii]
                fv[ii] = vg.fractg[ic2]
                vx[ii] = vg.ccell[ic2,1]; vy[ii] = vg.ccell[ic2,2]; vz[ii] = vg.ccell[ic2,3]
            end
            fv[nn+1] = vg.fractg[ic]
            vx[nn+1] = vg.ccell[ic,1]; vy[nn+1] = vg.ccell[ic,2]; vz[nn+1] = vg.ccell[ic,3]
            xn, yn, zn = lsrec(@view(fv[1:nn+1]), @view(vx[1:nn+1]),
                               @view(vy[1:nn+1]), @view(vz[1:nn+1]), pn)
        end
        vg.xnormg[ic] = xn; vg.ynormg[ic] = yn; vg.znormg[ic] = zn
        c = _enforce_volume(vg, ic, poly, xn, yn, zn; igrid)
        vg.rholig[ic] = c
    end
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════
#  2.  ELCIR – Extended Level-Contour Interface Reconstruction (irec=2)
# ═══════════════════════════════════════════════════════════════════════════
"""
    elcir!(vg::VOFGrid; iw=1, pn=1.5, igrid=1)

PLIC reconstruction using the extended level-contour (ELCIR).
"""
function elcir!(vg::VOFGrid; iw::Int=1, pn::Float64=1.5, igrid::Int=1)
    g      = vg.grid
    ncell  = g.ncell
    icint  = vg.icint
    work   = _get_recon_work(vg)
    iplc      = work.iplc
    nlc_v     = work.nlc_v
    ineigblc  = work.ineigblc
    xlc       = work.xlc
    ylc       = work.ylc
    zlc       = work.zlc
    iso_ws    = work.iso_ws
    iso_cells = work.iso_cells
    fill!(iplc, 0)
    fill!(nlc_v, 0)

    phivlc = 0.5
    # Local LC construction (over all cells)
    for ic in 1:ncell
        fnmax = 0.0; fnmin = 1.0
        nodes = g.ipcell[ic]
        nnodes = length(nodes)
        phi_v = @view(work.phi_v[1:nnodes])
        @inbounds for i in 1:nnodes
            phi_v[i] = vg.frnod[nodes[i]]
        end
        for p in phi_v; fnmax = max(fnmax, p); fnmin = min(fnmin, p); end
        (fnmax ≤ 0.5 || fnmin ≥ 0.5) && continue

        isoap!(iso_ws, iso_cells[ic], phi_v, phivlc)
        iso_ws.npoly == 1 || continue

        iplc[ic] = ic
        n_lc = iso_ws.poly_lens[1]
        nlc_v[ic] = n_lc
        neig_lc = _ensure_lc_storage!(ineigblc[ic], n_lc)
        xarr = _ensure_lc_storage!(xlc[ic], n_lc + 1)
        yarr = _ensure_lc_storage!(ylc[ic], n_lc + 1)
        zarr = _ensure_lc_storage!(zlc[ic], n_lc + 1)
        xsum = ysum = zsum = 0.0
        for i in 1:n_lc
            i2 = (i == n_lc) ? 1 : i + 1
            ip = iso_ws.polys[1, i]
            ip2 = iso_ws.polys[1, i2]
            is_face = iso_ws.isoeface[ip]
            iface = g.iscell[ic][is_face]
            icneig = g.icface[iface,1]
            icneig == ic && (icneig = g.icface[iface,2])
            neig_lc[i] = icneig
            mx = 0.5*(iso_ws.vertiso[ip][1] + iso_ws.vertiso[ip2][1])
            my = 0.5*(iso_ws.vertiso[ip][2] + iso_ws.vertiso[ip2][2])
            mz = 0.5*(iso_ws.vertiso[ip][3] + iso_ws.vertiso[ip2][3])
            xarr[i] = mx; yarr[i] = my; zarr[i] = mz
            xsum += mx; ysum += my; zsum += mz
        end
        xarr[end] = xsum/n_lc; yarr[end] = ysum/n_lc; zarr[end] = zsum/n_lc
    end

    # Preliminary normal + extended LC update
    for ind in eachindex(icint)
        ic = icint[ind]
        iplc[ic] == 0 && continue
        n_lc = nlc_v[ic]
        vertc = _load_vertc!(work.vertc, xlc[ic], ylc[ic], zlc[ic], n_lc + 1)
        xn, yn, zn = isrec(iw, vertc)
        vg.xnormg[ic] = xn; vg.ynormg[ic] = yn; vg.znormg[ic] = zn
        # Replace edge midpoints with neighbour LC centres
        for i in 1:n_lc
            icn = ineigblc[ic][i]
            (icn == 0 || iplc[icn] == 0) && continue
            nn = nlc_v[icn]
            xlc[ic][i] = xlc[icn][nn+1]
            ylc[ic][i] = ylc[icn][nn+1]
            zlc[ic][i] = zlc[icn][nn+1]
        end
    end

    # Volume enforcement
    for ind in eachindex(icint)
        ic = icint[ind]
        poly = defcell!(work.poly, work.ipg, vg, ic)

        if iplc[ic] != 0
            n_lc = nlc_v[ic]
            vertc = _load_vertc!(work.vertc, xlc[ic], ylc[ic], zlc[ic], n_lc + 1)
            xn, yn, zn = isrec(iw, vertc)
            α = acos(clamp(vg.xnormg[ic]*xn + vg.ynormg[ic]*yn + vg.znormg[ic]*zn, -1.0, 1.0))
            if α > 1.2
                xn, yn, zn = vg.xnormg[ic], vg.ynormg[ic], vg.znormg[ic]
            end
        else
            neigs = vg.ineigb[ic]; nn = length(neigs)
            fv = work.fv
            vx = work.vx
            vy = work.vy
            vz = work.vz
            for ii in 1:nn
                ic2 = neigs[ii]; fv[ii] = vg.fractg[ic2]
                vx[ii] = vg.ccell[ic2,1]; vy[ii] = vg.ccell[ic2,2]; vz[ii] = vg.ccell[ic2,3]
            end
            fv[nn+1] = vg.fractg[ic]; vx[nn+1] = vg.ccell[ic,1]; vy[nn+1] = vg.ccell[ic,2]; vz[nn+1] = vg.ccell[ic,3]
            xn, yn, zn = lsrec(@view(fv[1:nn+1]), @view(vx[1:nn+1]),
                               @view(vy[1:nn+1]), @view(vz[1:nn+1]), pn)
        end
        vg.xnormg[ic] = xn; vg.ynormg[ic] = yn; vg.znormg[ic] = zn
        c = _enforce_volume(vg, ic, poly, xn, yn, zn; igrid)
        vg.rholig[ic] = c
    end
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════
#  3.  CLCIR – Conservative Level-Contour Interface Reconstruction (irec=1)
# ═══════════════════════════════════════════════════════════════════════════
"""
    clcir!(vg::VOFGrid; iw=1, pn=1.5, igrid=1)

PLIC reconstruction using the conservative level-contour (CLCIR).
Same as ELCIR with a second-pass PLIC centroid update.
"""
function clcir!(vg::VOFGrid; iw::Int=1, pn::Float64=1.5, igrid::Int=1)
    # Perform ELCIR first
    elcir!(vg; iw, pn, igrid)

    g     = vg.grid
    icint = vg.icint
    work  = _get_recon_work(vg)

    # Compute PLIC centroids from current reconstruction
    cplic = work.cplic
    mark  = work.mark
    fill!(mark, 0)
    for ind in eachindex(icint)
        ic = icint[ind]
        poly = defcell!(work.poly, work.ipg, vg, ic)
        xnc, ync, znc = vg.xnormg[ic], vg.ynormg[ic], vg.znormg[ic]
        c = vg.rholig[ic]
        nts_before = poly.nts
        icontn, icontp = inte3d!(poly, c, xnc, ync, znc)
        ctr = _plic_centroid(poly, nts_before)
        if ctr !== nothing
            cplic[ic,1] = ctr[1]
            cplic[ic,2] = ctr[2]
            cplic[ic,3] = ctr[3]
            mark[ic] = 1
        end
    end

    # Second-pass: build extended LC using PLIC centroids and re-enforce
    # For cells with mark==1, replace neighbour LC centres with PLIC centres
    # then redo ISREC + enforce.  This is the "conservative" correction.
    for ind in eachindex(icint)
        ic = icint[ind]
        mark[ic] == 0 && continue
        poly = defcell!(work.poly, work.ipg, vg, ic)

        # Gather neighbour PLIC centroids
        neigs = vg.ineigb[ic]
        vertc = work.vertc
        empty!(vertc)
        for icn in neigs
            mark[icn] == 0 && continue
            push!(vertc, (cplic[icn,1], cplic[icn,2], cplic[icn,3]))
        end
        push!(vertc, (cplic[ic,1], cplic[ic,2], cplic[ic,3]))

        if length(vertc) >= 4
            xn, yn, zn = isrec(iw, vertc)
            α = acos(clamp(vg.xnormg[ic]*xn + vg.ynormg[ic]*yn + vg.znormg[ic]*zn, -1.0, 1.0))
            if α < 1.2
                vg.xnormg[ic] = xn; vg.ynormg[ic] = yn; vg.znormg[ic] = zn
                c2 = _enforce_volume(vg, ic, poly, xn, yn, zn; igrid)
                vg.rholig[ic] = c2
            end
        end
    end
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════
#  4.  LSGIR – Least-Squares Gradient Interface Reconstruction (irec=4)
# ═══════════════════════════════════════════════════════════════════════════
"""
    lsgir!(vg::VOFGrid; pn=1.5, igrid=1)

PLIC reconstruction using a weighted least-squares gradient estimate.
"""
function lsgir!(vg::VOFGrid; pn::Float64=1.5, igrid::Int=1)
    icint = vg.icint
    work = _get_recon_work(vg)
    for ind in eachindex(icint)
        ic = icint[ind]
        poly = defcell!(work.poly, work.ipg, vg, ic)
        neigs = vg.ineigb[ic]; nn = length(neigs)
        fv = work.fv
        vx = work.vx
        vy = work.vy
        vz = work.vz
        for ii in 1:nn
            ic2 = neigs[ii]; fv[ii] = vg.fractg[ic2]
            vx[ii] = vg.ccell[ic2,1]; vy[ii] = vg.ccell[ic2,2]; vz[ii] = vg.ccell[ic2,3]
        end
        fv[nn+1] = vg.fractg[ic]
        vx[nn+1] = vg.ccell[ic,1]; vy[nn+1] = vg.ccell[ic,2]; vz[nn+1] = vg.ccell[ic,3]
        xn, yn, zn = lsrec(@view(fv[1:nn+1]), @view(vx[1:nn+1]),
                           @view(vy[1:nn+1]), @view(vz[1:nn+1]), pn)
        vg.xnormg[ic] = xn; vg.ynormg[ic] = yn; vg.znormg[ic] = zn
        c = _enforce_volume(vg, ic, poly, xn, yn, zn; igrid)
        vg.rholig[ic] = c
    end
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════
#  5.  SWIR – Swartz Interface Reconstruction (irec=5)
# ═══════════════════════════════════════════════════════════════════════════
"""
    swrec(vg, ic, cplic, markplic; igrid) -> (xn, yn, zn)

Inner Swartz procedure for a single cell.
"""
function swrec(vg::VOFGrid, ic::Int, cplic::Matrix{Float64},
               markplic::Vector{Int}, work::_ReconWork; igrid::Int=1)
    neigs = vg.ineigb[ic]
    xnsw = ynsw = znsw = 0.0
    nsw = 0

    for ic2 in neigs
        markplic[ic2] == 0 && continue
        dp = vg.xnormg[ic]*vg.xnormg[ic2] + vg.ynormg[ic]*vg.ynormg[ic2] + vg.znormg[ic]*vg.znormg[ic2]
        α = acos(clamp(dp, -1.0, 1.0))
        α ≥ 0.785398 && continue  # > 45°

        xv = cplic[ic2,1] - cplic[ic,1]
        yv = cplic[ic2,2] - cplic[ic,2]
        zv = cplic[ic2,3] - cplic[ic,3]
        # Double cross-product
        nx0, ny0, nz0 = vg.xnormg[ic], vg.ynormg[ic], vg.znormg[ic]
        xv2 = yv*nz0 - zv*ny0; yv2 = zv*nx0 - xv*nz0; zv2 = xv*ny0 - yv*nx0
        xn = -yv*zv2 + zv*yv2; yn = -zv*xv2 + xv*zv2; zn = -xv*yv2 + yv*xv2

        s1 = xn*nx0 + yn*ny0 + zn*nz0
        s2 = xn*vg.xnormg[ic2] + yn*vg.ynormg[ic2] + zn*vg.znormg[ic2]
        d  = sqrt(xn^2 + yn^2 + zn^2)
        (s1 ≤ 0 || s2 ≤ 0 || d == 0) && continue

        xn /= d; yn /= d; zn /= d

        # Newton-like iterations on pairs of PLIC centroids
        for iter in 1:10
            xn0, yn0, zn0 = xn, yn, zn
            # Compute PLIC centroids for IC and IC2 with normal (xn0,yn0,zn0)
            poly1 = defcell!(work.poly1, work.ipg, vg, ic)
            c1 = _enforce_volume(vg, ic, poly1, xn0, yn0, zn0; igrid)
            nts1 = poly1.nts
            inte3d!(poly1, c1, xn0, yn0, zn0)
            ct1 = _plic_centroid(poly1, nts1)
            ct1 === nothing && (ct1 = (cplic[ic,1], cplic[ic,2], cplic[ic,3]))

            poly2 = defcell!(work.poly2, work.ipg, vg, ic2)
            v2 = vg.fractg[ic2] * vg.vcell[ic2]; vt2 = vg.vcell[ic2]
            if igrid == 1
                c2 = enforv3dsz(vg.boxcell[ic2,2]-vg.boxcell[ic2,1],
                                vg.boxcell[ic2,4]-vg.boxcell[ic2,3],
                                vg.boxcell[ic2,6]-vg.boxcell[ic2,5],
                                v2, poly2.vertp, xn0, yn0, zn0)
            else
                c2 = enforv3d(poly2, v2, vt2, xn0, yn0, zn0)
            end
            nts2 = poly2.nts
            inte3d!(poly2, c2, xn0, yn0, zn0)
            ct2 = _plic_centroid(poly2, nts2)
            ct2 === nothing && (ct2 = (cplic[ic2,1], cplic[ic2,2], cplic[ic2,3]))

            xvv = ct2[1]-ct1[1]; yvv = ct2[2]-ct1[2]; zvv = ct2[3]-ct1[3]
            xvv2 = yvv*zn0 - zvv*yn0; yvv2 = zvv*xn0 - xvv*zn0; zvv2 = xvv*yn0 - yvv*xn0
            xn = -yvv*zvv2 + zvv*yvv2; yn = -zvv*xvv2 + xvv*zvv2; zn = -xvv*yvv2 + yvv*xvv2
            dd = sqrt(xn^2 + yn^2 + zn^2)
            dd > 0 && (xn /= dd; yn /= dd; zn /= dd)
            toli = abs(acos(clamp(xn*xn0 + yn*yn0 + zn*zn0, -1.0, 1.0)))
            toli < 1e-6 && break
        end

        xnsw += xn; ynsw += yn; znsw += zn; nsw += 1
    end

    if nsw > 0
        d = sqrt(xnsw^2 + ynsw^2 + znsw^2)
        d > 0 && return (xnsw/d, ynsw/d, znsw/d)
    end
    return (vg.xnormg[ic], vg.ynormg[ic], vg.znormg[ic])
end

"""
    swir!(vg::VOFGrid; pn=1.5, niter=4, tolir=0.001, igrid=1)

PLIC reconstruction using the Swartz iterative algorithm.
"""
function swir!(vg::VOFGrid; pn::Float64=1.5, niter::Int=4, tolir::Float64=0.001, igrid::Int=1)
    icint = vg.icint
    work  = _get_recon_work(vg)

    cplic = work.cplic
    cplic2 = work.cplic2
    markplic = work.mark
    markplicend = work.mark2
    xnorm2 = work.xnorm2
    ynorm2 = work.ynorm2
    znorm2 = work.znorm2
    rholi2 = work.rholi2
    alpha = work.alpha
    fill!(markplic, 0)
    fill!(markplicend, 0)

    # Init with LSGIR
    for ind in eachindex(icint)
        ic = icint[ind]
        poly = defcell!(work.poly, work.ipg, vg, ic)
        neigs = vg.ineigb[ic]; nn = length(neigs)
        fv = work.fv
        vx = work.vx
        vy = work.vy
        vz = work.vz
        for ii in 1:nn
            ic2 = neigs[ii]; fv[ii] = vg.fractg[ic2]
            vx[ii] = vg.ccell[ic2,1]; vy[ii] = vg.ccell[ic2,2]; vz[ii] = vg.ccell[ic2,3]
        end
        fv[nn+1] = vg.fractg[ic]; vx[nn+1] = vg.ccell[ic,1]; vy[nn+1] = vg.ccell[ic,2]; vz[nn+1] = vg.ccell[ic,3]
        xn, yn, zn = lsrec(@view(fv[1:nn+1]), @view(vx[1:nn+1]),
                           @view(vy[1:nn+1]), @view(vz[1:nn+1]), pn)
        vg.xnormg[ic] = xn; vg.ynormg[ic] = yn; vg.znormg[ic] = zn
        c = _enforce_volume(vg, ic, poly, xn, yn, zn; igrid)
        vg.rholig[ic] = c

        nts_before = poly.nts
        inte3d!(poly, c, xn, yn, zn)
        ctr = _plic_centroid(poly, nts_before)
        if ctr !== nothing
            cplic[ic,1] = ctr[1]
            cplic[ic,2] = ctr[2]
            cplic[ic,3] = ctr[3]
            markplic[ic] = 1
        end
    end

    for iter in 1:niter
        alphalim = 0.523599 / iter  # 30°/iter

        for ind in eachindex(icint)
            ic = icint[ind]
            (markplic[ic] == 0 || markplicend[ic] != 0) && continue
            xn, yn, zn = swrec(vg, ic, cplic, markplic, work; igrid)
            xnorm2[ic] = xn; ynorm2[ic] = yn; znorm2[ic] = zn
            alpha[ic] = acos(clamp(vg.xnormg[ic]*xn + vg.ynormg[ic]*yn + vg.znormg[ic]*zn, -1.0, 1.0))
            alpha[ic] < tolir && (markplicend[ic] = 1)
            alpha[ic] < alphalim || continue
            poly = defcell!(work.poly, work.ipg, vg, ic)
            c2 = _enforce_volume(vg, ic, poly, xn, yn, zn; igrid)
            rholi2[ic] = c2
            if iter < niter
                nts_b = poly.nts
                inte3d!(poly, c2, xn, yn, zn)
                ctr = _plic_centroid(poly, nts_b)
                if ctr !== nothing
                    cplic2[ic,1] = ctr[1]
                    cplic2[ic,2] = ctr[2]
                    cplic2[ic,3] = ctr[3]
                else
                    cplic2[ic,1] = cplic[ic,1]
                    cplic2[ic,2] = cplic[ic,2]
                    cplic2[ic,3] = cplic[ic,3]
                end
            end
        end

        # Update
        for ind in eachindex(icint)
            ic = icint[ind]
            markplic[ic] == 0 && continue
            alpha[ic] < alphalim || continue
            vg.xnormg[ic] = xnorm2[ic]; vg.ynormg[ic] = ynorm2[ic]; vg.znormg[ic] = znorm2[ic]
            vg.rholig[ic] = rholi2[ic]
            if iter < niter
                cplic[ic,1] = cplic2[ic,1]
                cplic[ic,2] = cplic2[ic,2]
                cplic[ic,3] = cplic2[ic,3]
            end
        end
    end
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════
#  6.  LSFIR – Least-Squares Fit Interface Reconstruction (irec=6)
# ═══════════════════════════════════════════════════════════════════════════
"""
    lsfir!(vg::VOFGrid; pn=1.5, niter=4, tolir=0.001, igrid=1)

PLIC reconstruction using least-squares fitting of PLIC centroids.
"""
function lsfir!(vg::VOFGrid; pn::Float64=1.5, niter::Int=4, tolir::Float64=0.001, igrid::Int=1)
    icint = vg.icint
    work  = _get_recon_work(vg)

    cplic = work.cplic
    cplic2 = work.cplic2
    markplic = work.mark
    markplicend = work.mark2
    xnorm2 = work.xnorm2
    ynorm2 = work.ynorm2
    znorm2 = work.znorm2
    rholi2 = work.rholi2
    alpha = work.alpha
    fill!(markplic, 0)
    fill!(markplicend, 0)

    # Init with LSGIR
    for ind in eachindex(icint)
        ic = icint[ind]
        poly = defcell!(work.poly, work.ipg, vg, ic)
        neigs = vg.ineigb[ic]; nn = length(neigs)
        fv = work.fv
        vx = work.vx
        vy = work.vy
        vz = work.vz
        for ii in 1:nn
            ic2 = neigs[ii]; fv[ii] = vg.fractg[ic2]
            vx[ii] = vg.ccell[ic2,1]; vy[ii] = vg.ccell[ic2,2]; vz[ii] = vg.ccell[ic2,3]
        end
        fv[nn+1] = vg.fractg[ic]; vx[nn+1] = vg.ccell[ic,1]; vy[nn+1] = vg.ccell[ic,2]; vz[nn+1] = vg.ccell[ic,3]
        xn, yn, zn = lsrec(@view(fv[1:nn+1]), @view(vx[1:nn+1]),
                           @view(vy[1:nn+1]), @view(vz[1:nn+1]), pn)
        vg.xnormg[ic] = xn; vg.ynormg[ic] = yn; vg.znormg[ic] = zn
        c = _enforce_volume(vg, ic, poly, xn, yn, zn; igrid)
        vg.rholig[ic] = c

        nts_before = poly.nts
        inte3d!(poly, c, xn, yn, zn)
        ctr = _plic_centroid(poly, nts_before)
        if ctr !== nothing
            cplic[ic,1] = ctr[1]
            cplic[ic,2] = ctr[2]
            cplic[ic,3] = ctr[3]
            markplic[ic] = 1
        end
    end

    for iter in 1:niter
        alphalim = 0.523599 / iter

        for ind in eachindex(icint)
            ic = icint[ind]
            (markplic[ic] == 0 || markplicend[ic] != 0) && continue

            neigs = vg.ineigb[ic]
            n = 0
            vcx = work.vcx; vcy = work.vcy; vcz = work.vcz
            nnx = work.nnx; nny = work.nny; nnz = work.nnz
            for icn in neigs
                markplic[icn] == 1 || continue
                n += 1
                vcx[n] = cplic[icn,1]; vcy[n] = cplic[icn,2]; vcz[n] = cplic[icn,3]
                nnx[n] = vg.xnormg[icn]; nny[n] = vg.ynormg[icn]; nnz[n] = vg.znormg[icn]
            end
            vcx[n+1] = cplic[ic,1]; vcy[n+1] = cplic[ic,2]; vcz[n+1] = cplic[ic,3]
            nnx[n+1] = vg.xnormg[ic]; nny[n+1] = vg.ynormg[ic]; nnz[n+1] = vg.znormg[ic]

            xn, yn, zn = lsfrec(@view(vcx[1:n+1]), @view(vcy[1:n+1]), @view(vcz[1:n+1]),
                                 @view(nnx[1:n+1]), @view(nny[1:n+1]), @view(nnz[1:n+1]))
            xnorm2[ic] = xn; ynorm2[ic] = yn; znorm2[ic] = zn
            alpha[ic] = acos(clamp(vg.xnormg[ic]*xn + vg.ynormg[ic]*yn + vg.znormg[ic]*zn, -1.0, 1.0))
            alpha[ic] < tolir && (markplicend[ic] = 1)
            alpha[ic] < alphalim || continue
            poly = defcell!(work.poly, work.ipg, vg, ic)
            c2 = _enforce_volume(vg, ic, poly, xn, yn, zn; igrid)
            rholi2[ic] = c2
            if iter < niter
                nts_b = poly.nts
                inte3d!(poly, c2, xn, yn, zn)
                ctr = _plic_centroid(poly, nts_b)
                if ctr !== nothing
                    cplic2[ic,1] = ctr[1]
                    cplic2[ic,2] = ctr[2]
                    cplic2[ic,3] = ctr[3]
                else
                    cplic2[ic,1] = cplic[ic,1]
                    cplic2[ic,2] = cplic[ic,2]
                    cplic2[ic,3] = cplic[ic,3]
                end
            end
        end

        for ind in eachindex(icint)
            ic = icint[ind]
            markplic[ic] == 0 && continue
            alpha[ic] < alphalim || continue
            vg.xnormg[ic] = xnorm2[ic]; vg.ynormg[ic] = ynorm2[ic]; vg.znormg[ic] = znorm2[ic]
            vg.rholig[ic] = rholi2[ic]
            if iter < niter
                cplic[ic,1] = cplic2[ic,1]
                cplic[ic,2] = cplic2[ic,2]
                cplic[ic,3] = cplic2[ic,3]
            end
        end
    end
    return nothing
end
