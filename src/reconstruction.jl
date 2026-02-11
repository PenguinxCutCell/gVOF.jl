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
    # Build weighted AᵀA and Aᵀb
    ata = zeros(3, 3)
    atb = zeros(3)
    xr, yr, zr, fr = vertcx[n+1], vertcy[n+1], vertcz[n+1], f[n+1]
    for i in 1:n
        d = sqrt((vertcx[i]-xr)^2 + (vertcy[i]-yr)^2 + (vertcz[i]-zr)^2)
        w = 1.0 / d^pn
        a1 = w * (vertcx[i] - xr)
        a2 = w * (vertcy[i] - yr)
        a3 = w * (vertcz[i] - zr)
        b  = w * (f[i] - fr)
        a  = (a1, a2, a3)
        for j in 1:3, k in 1:3
            ata[j,k] += a[j] * a[k]
        end
        for j in 1:3
            atb[j] += a[j] * b
        end
    end

    # Solve 3×3 system (Cramer / back-substitution with pivoting guard)
    xn = yn = zn = 0.0
    try
        sol = ata \ atb
        xn, yn, zn = sol[1], sol[2], sol[3]
    catch
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

# ═══════════════════════════════════════════════════════════════════════════
#  1.  LLCIR – Local Level-Contour Interface Reconstruction (irec=3)
# ═══════════════════════════════════════════════════════════════════════════
"""
    llcir!(vg::VOFGrid; iw=1, pn=1.5, igrid=1)

PLIC reconstruction using the local level-contour iso-surface (LLCIR).
"""
function llcir!(vg::VOFGrid; iw::Int=1, pn::Float64=1.5, igrid::Int=1)
    g      = vg.grid
    ncell  = g.ncell
    icint  = vg.icint

    # Per-cell LC data
    iplc  = zeros(Int, ncell)
    nlc_v = zeros(Int, ncell)
    xlc   = Vector{Vector{Float64}}(undef, ncell)
    ylc   = Vector{Vector{Float64}}(undef, ncell)
    zlc   = Vector{Vector{Float64}}(undef, ncell)

    # Local iso-surface construction
    phivlc = 0.5
    for ind in eachindex(icint)
        ic = icint[ind]
        fnmax = 0.0; fnmin = 1.0
        nodes = g.ipcell[ic]
        phi_v = Float64[vg.frnod[ip] for ip in nodes]
        for p in phi_v
            fnmax = max(fnmax, p); fnmin = min(fnmin, p)
        end
        fnmax ≤ 0.5 && continue
        fnmin ≥ 0.5 && continue

        poly_isoap = cellgrid(g, ic)
        res = isoap(poly_isoap, phi_v, phivlc)
        niso(res) == 1 || continue

        iplc[ic] = ic
        ipv_iso = res.ipviso[1]
        n_lc = length(ipv_iso)
        nlc_v[ic] = n_lc
        xarr = zeros(n_lc + 1); yarr = zeros(n_lc + 1); zarr = zeros(n_lc + 1)
        xsum = ysum = zsum = 0.0
        for i in 1:n_lc
            ip  = ipv_iso[i]
            v   = res.vertiso[ip]
            xarr[i] = v[1]; yarr[i] = v[2]; zarr[i] = v[3]
            xsum += v[1]; ysum += v[2]; zsum += v[3]
        end
        xarr[end] = xsum / n_lc; yarr[end] = ysum / n_lc; zarr[end] = zsum / n_lc
        xlc[ic] = xarr; ylc[ic] = yarr; zlc[ic] = zarr
    end

    # Volume enforcement
    for ind in eachindex(icint)
        ic = icint[ind]
        poly = defcell(vg, ic)

        if iplc[ic] != 0
            n_lc = nlc_v[ic]
            vertc = NTuple{3,Float64}[(xlc[ic][i], ylc[ic][i], zlc[ic][i]) for i in 1:n_lc+1]
            xn, yn, zn = isrec(iw, vertc)
        else
            neigs = vg.ineigb[ic]
            nn    = length(neigs)
            fv    = zeros(nn+1); vx = zeros(nn+1); vy = zeros(nn+1); vz = zeros(nn+1)
            for ii in 1:nn
                ic2 = neigs[ii]
                fv[ii] = vg.fractg[ic2]
                vx[ii] = vg.ccell[ic2,1]; vy[ii] = vg.ccell[ic2,2]; vz[ii] = vg.ccell[ic2,3]
            end
            fv[nn+1] = vg.fractg[ic]
            vx[nn+1] = vg.ccell[ic,1]; vy[nn+1] = vg.ccell[ic,2]; vz[nn+1] = vg.ccell[ic,3]
            xn, yn, zn = lsrec(fv, vx, vy, vz, pn)
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

    iplc      = zeros(Int, ncell)
    nlc_v     = zeros(Int, ncell)
    ineigblc  = Vector{Vector{Int}}(undef, ncell)
    xlc = Vector{Vector{Float64}}(undef, ncell)
    ylc = Vector{Vector{Float64}}(undef, ncell)
    zlc = Vector{Vector{Float64}}(undef, ncell)

    phivlc = 0.5
    # Local LC construction (over all cells)
    for ic in 1:ncell
        iplc[ic] = 0
        fnmax = 0.0; fnmin = 1.0
        nodes = g.ipcell[ic]
        phi_v = Float64[vg.frnod[ip] for ip in nodes]
        for p in phi_v; fnmax = max(fnmax, p); fnmin = min(fnmin, p); end
        (fnmax ≤ 0.5 || fnmin ≥ 0.5) && continue

        poly_is = cellgrid(g, ic)
        res = isoap(poly_is, phi_v, phivlc)
        niso(res) == 1 || continue

        iplc[ic] = ic
        ipv_iso = res.ipviso[1]
        n_lc = length(ipv_iso)
        nlc_v[ic] = n_lc
        neig_lc = Int[]
        xarr = zeros(n_lc + 1); yarr = zeros(n_lc + 1); zarr = zeros(n_lc + 1)
        xsum = ysum = zsum = 0.0
        for i in 1:n_lc
            i2 = (i == n_lc) ? 1 : i + 1
            ip = ipv_iso[i]; ip2 = ipv_iso[i2]
            is_face = res.isoeface[ip]
            iface = g.iscell[ic][is_face]
            icneig = g.icface[iface,1]
            icneig == ic && (icneig = g.icface[iface,2])
            push!(neig_lc, icneig)
            mx = 0.5*(res.vertiso[ip][1] + res.vertiso[ip2][1])
            my = 0.5*(res.vertiso[ip][2] + res.vertiso[ip2][2])
            mz = 0.5*(res.vertiso[ip][3] + res.vertiso[ip2][3])
            xarr[i] = mx; yarr[i] = my; zarr[i] = mz
            xsum += mx; ysum += my; zsum += mz
        end
        xarr[end] = xsum/n_lc; yarr[end] = ysum/n_lc; zarr[end] = zsum/n_lc
        ineigblc[ic] = neig_lc
        xlc[ic] = xarr; ylc[ic] = yarr; zlc[ic] = zarr
    end

    # Preliminary normal + extended LC update
    for ind in eachindex(icint)
        ic = icint[ind]
        iplc[ic] == 0 && continue
        n_lc = nlc_v[ic]
        vertc = NTuple{3,Float64}[(xlc[ic][i], ylc[ic][i], zlc[ic][i]) for i in 1:n_lc+1]
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
        poly = defcell(vg, ic)

        if iplc[ic] != 0
            n_lc = nlc_v[ic]
            vertc = NTuple{3,Float64}[(xlc[ic][i], ylc[ic][i], zlc[ic][i]) for i in 1:n_lc+1]
            xn, yn, zn = isrec(iw, vertc)
            α = acos(clamp(vg.xnormg[ic]*xn + vg.ynormg[ic]*yn + vg.znormg[ic]*zn, -1.0, 1.0))
            if α > 1.2
                xn, yn, zn = vg.xnormg[ic], vg.ynormg[ic], vg.znormg[ic]
            end
        else
            neigs = vg.ineigb[ic]; nn = length(neigs)
            fv = zeros(nn+1); vx = zeros(nn+1); vy = zeros(nn+1); vz = zeros(nn+1)
            for ii in 1:nn
                ic2 = neigs[ii]; fv[ii] = vg.fractg[ic2]
                vx[ii] = vg.ccell[ic2,1]; vy[ii] = vg.ccell[ic2,2]; vz[ii] = vg.ccell[ic2,3]
            end
            fv[nn+1] = vg.fractg[ic]; vx[nn+1] = vg.ccell[ic,1]; vy[nn+1] = vg.ccell[ic,2]; vz[nn+1] = vg.ccell[ic,3]
            xn, yn, zn = lsrec(fv, vx, vy, vz, pn)
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
    ncell = g.ncell
    icint = vg.icint

    # Compute PLIC centroids from current reconstruction
    cplic = zeros(ncell, 3)
    mark  = zeros(Int, ncell)
    for ind in eachindex(icint)
        ic = icint[ind]
        poly = defcell(vg, ic)
        xnc, ync, znc = vg.xnormg[ic], vg.ynormg[ic], vg.znormg[ic]
        c = vg.rholig[ic]
        nts_before = poly.nts
        icontn, icontp = inte3d!(poly, c, xnc, ync, znc)
        ctr = _plic_centroid(poly, nts_before)
        if ctr !== nothing
            cplic[ic,:] .= ctr
            mark[ic] = 1
        end
    end

    # Second-pass: build extended LC using PLIC centroids and re-enforce
    # For cells with mark==1, replace neighbour LC centres with PLIC centres
    # then redo ISREC + enforce.  This is the "conservative" correction.
    for ind in eachindex(icint)
        ic = icint[ind]
        mark[ic] == 0 && continue
        poly = defcell(vg, ic)

        # Gather neighbour PLIC centroids
        neigs = vg.ineigb[ic]
        vertc = NTuple{3,Float64}[]
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
    for ind in eachindex(icint)
        ic = icint[ind]
        poly = defcell(vg, ic)
        neigs = vg.ineigb[ic]; nn = length(neigs)
        fv = zeros(nn+1); vx = zeros(nn+1); vy = zeros(nn+1); vz = zeros(nn+1)
        for ii in 1:nn
            ic2 = neigs[ii]; fv[ii] = vg.fractg[ic2]
            vx[ii] = vg.ccell[ic2,1]; vy[ii] = vg.ccell[ic2,2]; vz[ii] = vg.ccell[ic2,3]
        end
        fv[nn+1] = vg.fractg[ic]
        vx[nn+1] = vg.ccell[ic,1]; vy[nn+1] = vg.ccell[ic,2]; vz[nn+1] = vg.ccell[ic,3]
        xn, yn, zn = lsrec(fv, vx, vy, vz, pn)
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
               markplic::Vector{Int}; igrid::Int=1)
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
            poly1 = defcell(vg, ic)
            c1 = _enforce_volume(vg, ic, poly1, xn0, yn0, zn0; igrid)
            nts1 = poly1.nts
            inte3d!(poly1, c1, xn0, yn0, zn0)
            ct1 = _plic_centroid(poly1, nts1)
            ct1 === nothing && (ct1 = (cplic[ic,1], cplic[ic,2], cplic[ic,3]))

            poly2 = defcell(vg, ic2)
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
    g     = vg.grid
    ncell = g.ncell
    icint = vg.icint

    cplic = zeros(ncell, 3)
    markplic    = zeros(Int, ncell)
    markplicend = zeros(Int, ncell)

    # Init with LSGIR
    for ind in eachindex(icint)
        ic = icint[ind]
        poly = defcell(vg, ic)
        neigs = vg.ineigb[ic]; nn = length(neigs)
        fv = zeros(nn+1); vx = zeros(nn+1); vy = zeros(nn+1); vz = zeros(nn+1)
        for ii in 1:nn
            ic2 = neigs[ii]; fv[ii] = vg.fractg[ic2]
            vx[ii] = vg.ccell[ic2,1]; vy[ii] = vg.ccell[ic2,2]; vz[ii] = vg.ccell[ic2,3]
        end
        fv[nn+1] = vg.fractg[ic]; vx[nn+1] = vg.ccell[ic,1]; vy[nn+1] = vg.ccell[ic,2]; vz[nn+1] = vg.ccell[ic,3]
        xn, yn, zn = lsrec(fv, vx, vy, vz, pn)
        vg.xnormg[ic] = xn; vg.ynormg[ic] = yn; vg.znormg[ic] = zn
        c = _enforce_volume(vg, ic, poly, xn, yn, zn; igrid)
        vg.rholig[ic] = c

        nts_before = poly.nts
        inte3d!(poly, c, xn, yn, zn)
        ctr = _plic_centroid(poly, nts_before)
        if ctr !== nothing
            cplic[ic,:] .= ctr; markplic[ic] = 1
        end
    end

    # Swartz iterations
    xnorm2 = zeros(ncell); ynorm2 = zeros(ncell); znorm2 = zeros(ncell)
    rholi2 = zeros(ncell); alpha  = zeros(ncell)
    cplic2 = zeros(ncell, 3)

    for iter in 1:niter
        alphalim = 0.523599 / iter  # 30°/iter

        for ind in eachindex(icint)
            ic = icint[ind]
            (markplic[ic] == 0 || markplicend[ic] != 0) && continue
            xn, yn, zn = swrec(vg, ic, cplic, markplic; igrid)
            xnorm2[ic] = xn; ynorm2[ic] = yn; znorm2[ic] = zn
            alpha[ic] = acos(clamp(vg.xnormg[ic]*xn + vg.ynormg[ic]*yn + vg.znormg[ic]*zn, -1.0, 1.0))
            alpha[ic] < tolir && (markplicend[ic] = 1)
            alpha[ic] < alphalim || continue
            poly = defcell(vg, ic)
            c2 = _enforce_volume(vg, ic, poly, xn, yn, zn; igrid)
            rholi2[ic] = c2
            if iter < niter
                nts_b = poly.nts
                inte3d!(poly, c2, xn, yn, zn)
                ctr = _plic_centroid(poly, nts_b)
                cplic2[ic,:] .= ctr !== nothing ? ctr : cplic[ic,:]
            end
        end

        # Update
        for ind in eachindex(icint)
            ic = icint[ind]
            markplic[ic] == 0 && continue
            alpha[ic] < alphalim || continue
            vg.xnormg[ic] = xnorm2[ic]; vg.ynormg[ic] = ynorm2[ic]; vg.znormg[ic] = znorm2[ic]
            vg.rholig[ic] = rholi2[ic]
            iter < niter && (cplic[ic,:] .= cplic2[ic,:])
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
    g     = vg.grid
    ncell = g.ncell
    icint = vg.icint

    cplic = zeros(ncell, 3)
    markplic    = zeros(Int, ncell)
    markplicend = zeros(Int, ncell)

    # Init with LSGIR
    for ind in eachindex(icint)
        ic = icint[ind]
        poly = defcell(vg, ic)
        neigs = vg.ineigb[ic]; nn = length(neigs)
        fv = zeros(nn+1); vx = zeros(nn+1); vy = zeros(nn+1); vz = zeros(nn+1)
        for ii in 1:nn
            ic2 = neigs[ii]; fv[ii] = vg.fractg[ic2]
            vx[ii] = vg.ccell[ic2,1]; vy[ii] = vg.ccell[ic2,2]; vz[ii] = vg.ccell[ic2,3]
        end
        fv[nn+1] = vg.fractg[ic]; vx[nn+1] = vg.ccell[ic,1]; vy[nn+1] = vg.ccell[ic,2]; vz[nn+1] = vg.ccell[ic,3]
        xn, yn, zn = lsrec(fv, vx, vy, vz, pn)
        vg.xnormg[ic] = xn; vg.ynormg[ic] = yn; vg.znormg[ic] = zn
        c = _enforce_volume(vg, ic, poly, xn, yn, zn; igrid)
        vg.rholig[ic] = c

        nts_before = poly.nts
        inte3d!(poly, c, xn, yn, zn)
        ctr = _plic_centroid(poly, nts_before)
        if ctr !== nothing
            cplic[ic,:] .= ctr; markplic[ic] = 1
        end
    end

    xnorm2 = zeros(ncell); ynorm2 = zeros(ncell); znorm2 = zeros(ncell)
    rholi2 = zeros(ncell); alpha  = zeros(ncell)
    cplic2 = zeros(ncell, 3)

    for iter in 1:niter
        alphalim = 0.523599 / iter

        for ind in eachindex(icint)
            ic = icint[ind]
            (markplic[ic] == 0 || markplicend[ic] != 0) && continue

            neigs = vg.ineigb[ic]
            n = 0
            max_n = length(neigs)
            vcx = zeros(max_n+1); vcy = zeros(max_n+1); vcz = zeros(max_n+1)
            nnx = zeros(max_n+1); nny = zeros(max_n+1); nnz = zeros(max_n+1)
            for icn in neigs
                markplic[icn] == 1 || continue
                n += 1
                vcx[n] = cplic[icn,1]; vcy[n] = cplic[icn,2]; vcz[n] = cplic[icn,3]
                nnx[n] = vg.xnormg[icn]; nny[n] = vg.ynormg[icn]; nnz[n] = vg.znormg[icn]
            end
            vcx[n+1] = cplic[ic,1]; vcy[n+1] = cplic[ic,2]; vcz[n+1] = cplic[ic,3]
            nnx[n+1] = vg.xnormg[ic]; nny[n+1] = vg.ynormg[ic]; nnz[n+1] = vg.znormg[ic]

            xn, yn, zn = lsfrec(vcx[1:n+1], vcy[1:n+1], vcz[1:n+1],
                                 nnx[1:n+1], nny[1:n+1], nnz[1:n+1])
            xnorm2[ic] = xn; ynorm2[ic] = yn; znorm2[ic] = zn
            alpha[ic] = acos(clamp(vg.xnormg[ic]*xn + vg.ynormg[ic]*yn + vg.znormg[ic]*zn, -1.0, 1.0))
            alpha[ic] < tolir && (markplicend[ic] = 1)
            alpha[ic] < alphalim || continue
            poly = defcell(vg, ic)
            c2 = _enforce_volume(vg, ic, poly, xn, yn, zn; igrid)
            rholi2[ic] = c2
            if iter < niter
                nts_b = poly.nts
                inte3d!(poly, c2, xn, yn, zn)
                ctr = _plic_centroid(poly, nts_b)
                cplic2[ic,:] .= ctr !== nothing ? ctr : cplic[ic,:]
            end
        end

        for ind in eachindex(icint)
            ic = icint[ind]
            markplic[ic] == 0 && continue
            alpha[ic] < alphalim || continue
            vg.xnormg[ic] = xnorm2[ic]; vg.ynormg[ic] = ynorm2[ic]; vg.znormg[ic] = znorm2[ic]
            vg.rholig[ic] = rholi2[ic]
            iter < niter && (cplic[ic,:] .= cplic2[ic,:])
        end
    end
    return nothing
end
