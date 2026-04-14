# ---------------------------------------------------------------------------
#  3D equality-constrained LVIRA reconstruction
# ---------------------------------------------------------------------------

const _LVIRA3D_WEIGHT_P = 2.0

@inline function _canonicalize_normal(nx::Float64, ny::Float64, nz::Float64)
    d = sqrt(nx * nx + ny * ny + nz * nz)
    d > 0.0 || return (1.0, 0.0, 0.0)
    nx /= d
    ny /= d
    nz /= d

    eps = 1e-14
    if nz < -eps ||
       (abs(nz) <= eps && ny < -eps) ||
       (abs(nz) <= eps && abs(ny) <= eps && nx < 0.0)
        nx = -nx
        ny = -ny
        nz = -nz
    end
    return (nx, ny, nz)
end

@inline function _angles_to_normal(theta::Float64, phi::Float64)
    st = sin(theta)
    ct = cos(theta)
    sp = sin(phi)
    cp = cos(phi)
    nx = ct * sp
    ny = st * sp
    nz = cp
    d = sqrt(nx * nx + ny * ny + nz * nz)
    d > 0.0 || return (1.0, 0.0, 0.0)
    return (nx / d, ny / d, nz / d)
end

@inline function _normal_to_angles(nx::Float64, ny::Float64, nz::Float64)
    d = sqrt(nx * nx + ny * ny + nz * nz)
    d > 0.0 || return (0.0, π / 2.0)
    nx /= d
    ny /= d
    nz /= d
    theta = mod(atan(ny, nx), 2.0 * π)
    phi = acos(clamp(nz, -1.0, 1.0))
    phi = clamp(phi, 0.0, π)
    return (theta, phi)
end

@inline _wrap_theta(theta::Float64) = mod(theta, 2.0 * π)
@inline _clamp_phi(phi::Float64) = clamp(phi, 0.0, π)

@inline function _enforce_center_volume!(vg::VOFGrid, ic::Int, work::_ReconWork,
                                         nx::Float64, ny::Float64, nz::Float64;
                                         igrid::Int=1)
    poly = defcell!(work.poly, work.ipg, vg, ic)
    # Equality constraint: center-cell fraction must be matched exactly.
    d = _enforce_volume(vg, ic, poly, nx, ny, nz; igrid=igrid)
    return d
end

function _neighbor_stencil_ids!(vg::VOFGrid, ic::Int, work::_ReconWork;
                                mode::Symbol=:all_immediate,
                                tolfr::Float64=1e-12)
    mode == :all_immediate || error("Unsupported LVIRA3D stencil mode: $mode")

    ids = work.neigh_ids
    touched = work.neigh_touched
    mark = work.mark
    n = 0
    nt = 0

    for icn in vg.ineigb[ic]
        (icn == ic || icn < 1 || icn > vg.grid.ncell) && continue
        mark[icn] == 1 && continue

        if mode == :mixed_only
            f = vg.fractg[icn]
            (f <= tolfr || f >= 1.0 - tolfr) && continue
        end

        n += 1
        ids[n] = icn
        mark[icn] = 1
        nt += 1
        touched[nt] = icn
    end

    for k in 1:nt
        mark[touched[k]] = 0
    end

    return n
end

function _fill_neighbor_weights!(vg::VOFGrid, ic::Int, n_neigh::Int,
                                 work::_ReconWork; wtype::Symbol=:uniform)
    w = work.neigh_weights
    ids = work.neigh_ids
    if wtype == :uniform
        @inbounds for k in 1:n_neigh
            w[k] = 1.0
        end
    elseif wtype == :distance
        cx = vg.ccell[ic,1]
        cy = vg.ccell[ic,2]
        cz = vg.ccell[ic,3]
        @inbounds for k in 1:n_neigh
            icn = ids[k]
            dx = vg.ccell[icn,1] - cx
            dy = vg.ccell[icn,2] - cy
            dz = vg.ccell[icn,3] - cz
            d = max(1e-14, sqrt(dx * dx + dy * dy + dz * dz))
            w[k] = 1.0 / d^_LVIRA3D_WEIGHT_P
        end
    else
        error("Unsupported LVIRA3D weight type: $wtype")
    end
    return nothing
end

function _lsgir3d_cell!(vg::VOFGrid, ic::Int, work::_ReconWork;
                        pn::Float64=1.5, igrid::Int=1)
    neigs = vg.ineigb[ic]
    nn = length(neigs)
    if nn == 0
        nx, ny, nz = 1.0, 0.0, 0.0
        d = _enforce_center_volume!(vg, ic, work, nx, ny, nz; igrid=igrid)
        return nx, ny, nz, d
    end

    fv = work.fv
    vx = work.vx
    vy = work.vy
    vz = work.vz
    @inbounds for ii in 1:nn
        ic2 = neigs[ii]
        fv[ii] = vg.fractg[ic2]
        vx[ii] = vg.ccell[ic2,1]
        vy[ii] = vg.ccell[ic2,2]
        vz[ii] = vg.ccell[ic2,3]
    end
    fv[nn+1] = vg.fractg[ic]
    vx[nn+1] = vg.ccell[ic,1]
    vy[nn+1] = vg.ccell[ic,2]
    vz[nn+1] = vg.ccell[ic,3]

    nx, ny, nz = lsrec(@view(fv[1:nn+1]), @view(vx[1:nn+1]),
                       @view(vy[1:nn+1]), @view(vz[1:nn+1]), pn)
    d = _enforce_center_volume!(vg, ic, work, nx, ny, nz; igrid=igrid)
    return nx, ny, nz, d
end

function _lvira3d_seed_normal(vg::VOFGrid, ic::Int, work::_ReconWork;
                              seed::Symbol=:lsgir, pn::Float64=1.5)
    if seed == :current || seed == :swir
        nx = vg.xnormg[ic]
        ny = vg.ynormg[ic]
        nz = vg.znormg[ic]
        if isfinite(nx) && isfinite(ny) && isfinite(nz)
            d = sqrt(nx * nx + ny * ny + nz * nz)
            if d > 1e-14
                return (nx / d, ny / d, nz / d)
            end
        end
        # Fall through to LSGIR seed.
    elseif seed != :lsgir
        error("Unsupported LVIRA3D seed strategy: $seed")
    end

    neigs = vg.ineigb[ic]
    nn = length(neigs)
    if nn == 0
        return (1.0, 0.0, 0.0)
    end

    fv = work.fv
    vx = work.vx
    vy = work.vy
    vz = work.vz
    @inbounds for ii in 1:nn
        ic2 = neigs[ii]
        fv[ii] = vg.fractg[ic2]
        vx[ii] = vg.ccell[ic2,1]
        vy[ii] = vg.ccell[ic2,2]
        vz[ii] = vg.ccell[ic2,3]
    end
    fv[nn+1] = vg.fractg[ic]
    vx[nn+1] = vg.ccell[ic,1]
    vy[nn+1] = vg.ccell[ic,2]
    vz[nn+1] = vg.ccell[ic,3]

    nx, ny, nz = lsrec(@view(fv[1:nn+1]), @view(vx[1:nn+1]),
                       @view(vy[1:nn+1]), @view(vz[1:nn+1]), pn)
    d = sqrt(nx * nx + ny * ny + nz * nz)
    d > 0.0 || return (1.0, 0.0, 0.0)
    return (nx / d, ny / d, nz / d)
end

function _lvira3d_objective!(vg::VOFGrid, ic::Int, n_neigh::Int,
                             nx::Float64, ny::Float64, nz::Float64,
                             work::_ReconWork;
                             igrid::Int=1,
                             wtype::Symbol=:uniform)
    d = _enforce_center_volume!(vg, ic, work, nx, ny, nz; igrid=igrid)
    if !isfinite(d)
        return (Inf, d)
    end

    _fill_neighbor_weights!(vg, ic, n_neigh, work; wtype=wtype)

    ids = work.neigh_ids
    ws = work.neigh_weights
    j = 0.0
    @inbounds for k in 1:n_neigh
        icn = ids[k]
        # SAME plane (n,d) is extended to all stencil neighbors.
        fhat = _fraction_for_cell_with_plane!(work, vg, icn, nx, ny, nz, d)
        if !isfinite(fhat)
            return (Inf, d)
        end
        df = vg.fractg[icn] - fhat
        j += ws[k] * df * df
    end

    return (j, d)
end

function _lvira3d_optimize!(vg::VOFGrid, ic::Int, n_neigh::Int,
                            n0::NTuple{3,Float64}, work::_ReconWork;
                            igrid::Int=1,
                            tolir::Float64=1e-8,
                            niter::Int=80,
                            wtype::Symbol=:uniform)
    nx0, ny0, nz0 = n0
    d0 = sqrt(nx0 * nx0 + ny0 * ny0 + nz0 * nz0)
    if d0 <= 0.0
        nx0, ny0, nz0 = 1.0, 0.0, 0.0
    else
        nx0 /= d0
        ny0 /= d0
        nz0 /= d0
    end
    theta0, phi0 = _normal_to_angles(nx0, ny0, nz0)

    dtheta = 25.0 * π / 180.0
    dphi = 20.0 * π / 180.0

    ntheta = 11
    nphi = 9
    best_theta = theta0
    best_phi = phi0
    best_j = Inf
    best_d = NaN

    theta_min = theta0 - dtheta
    theta_max = theta0 + dtheta
    phi_min = max(0.0, phi0 - dphi)
    phi_max = min(π, phi0 + dphi)

    if phi_max <= phi_min
        phi_min = 0.0
        phi_max = π
    end

    for it in 0:(ntheta - 1)
        theta = theta_min + (theta_max - theta_min) * it / max(1, ntheta - 1)
        theta = _wrap_theta(theta)
        for ip in 0:(nphi - 1)
            phi = phi_min + (phi_max - phi_min) * ip / max(1, nphi - 1)
            phi = _clamp_phi(phi)
            nx, ny, nz = _angles_to_normal(theta, phi)
            j, d = _lvira3d_objective!(vg, ic, n_neigh, nx, ny, nz, work;
                                       igrid=igrid, wtype=wtype)
            if isfinite(j) && j < best_j
                best_j = j
                best_d = d
                best_theta = theta
                best_phi = phi
            end
        end
    end

    if !isfinite(best_j)
        return (nx0, ny0, nz0, NaN, Inf)
    end

    # Global fallback scan on the full sphere for robustness.
    # This is especially important when the local seed is poor.
    if best_j > 1e-6
        ntheta_g = 28
        nphi_g = 14
        for it in 0:(ntheta_g - 1)
            theta = 2.0 * π * it / max(1, ntheta_g - 1)
            for ip in 0:(nphi_g - 1)
                phi = π * ip / max(1, nphi_g - 1)
                nx, ny, nz = _angles_to_normal(theta, phi)
                j, d = _lvira3d_objective!(vg, ic, n_neigh, nx, ny, nz, work;
                                           igrid=igrid, wtype=wtype)
                if isfinite(j) && j < best_j
                    best_j = j
                    best_d = d
                    best_theta = theta
                    best_phi = phi
                end
            end
        end
    end

    htheta = 10.0 * π / 180.0
    hphi = 8.0 * π / 180.0

    for _ in 1:niter
        improved = false
        trial_pts = (
            (best_theta + htheta, best_phi),
            (best_theta - htheta, best_phi),
            (best_theta, best_phi + hphi),
            (best_theta, best_phi - hphi),
            (best_theta + htheta, best_phi + hphi),
            (best_theta + htheta, best_phi - hphi),
            (best_theta - htheta, best_phi + hphi),
            (best_theta - htheta, best_phi - hphi),
        )

        for pt in trial_pts
            theta = _wrap_theta(pt[1])
            phi = _clamp_phi(pt[2])
            nx, ny, nz = _angles_to_normal(theta, phi)
            j, d = _lvira3d_objective!(vg, ic, n_neigh, nx, ny, nz, work;
                                       igrid=igrid, wtype=wtype)
            if isfinite(j) && j < best_j - 1e-16
                best_j = j
                best_d = d
                best_theta = theta
                best_phi = phi
                improved = true
            end
        end

        if !improved
            htheta *= 0.5
            hphi *= 0.5
        end

        if max(htheta, hphi) < tolir
            break
        end
    end

    nx, ny, nz = _angles_to_normal(best_theta, best_phi)
    if !isfinite(best_d)
        best_d = _enforce_center_volume!(vg, ic, work, nx, ny, nz; igrid=igrid)
    end
    return (nx, ny, nz, best_d, best_j)
end

"""
    _elvira3d_candidates_placeholder(...)

Placeholder for a future optional canonical 3D ELVIRA implementation
(5×5×5 stencil with large candidate families on Cartesian grids).
"""
function _elvira3d_candidates_placeholder(args...)
    _ = args
    return nothing
end

"""
    lvira3d!(vg::VOFGrid; igrid=1, tolir=1e-10, niter=80, seed=:lsgir, wtype=:uniform)

3D equality-constrained LVIRA reconstruction for arbitrary polyhedral grids.
The method minimizes neighbor-cell mismatch while enforcing the center-cell
volume fraction exactly, using the same plane in all neighbors.

This method is generally more expensive than LSGIR but provides improved local
geometric consistency.
"""
function lvira3d!(vg::VOFGrid;
                  igrid::Int=1,
                  tolir::Float64=1e-10,
                  niter::Int=80,
                  seed::Symbol=:lsgir,
                  wtype::Symbol=:uniform)
    work = _get_recon_work(vg)
    tolfr = 1e-12

    for ic in vg.icint
        f = vg.fractg[ic]
        (f <= tolfr || f >= 1.0 - tolfr) && continue

        n_neigh = _neighbor_stencil_ids!(vg, ic, work; mode=:all_immediate, tolfr=tolfr)
        if n_neigh < 3
            nx, ny, nz, d = _lsgir3d_cell!(vg, ic, work; pn=1.5, igrid=igrid)
            vg.xnormg[ic] = nx
            vg.ynormg[ic] = ny
            vg.znormg[ic] = nz
            vg.rholig[ic] = d
            continue
        end

        n0 = _lvira3d_seed_normal(vg, ic, work; seed=seed, pn=1.5)
        nx, ny, nz, d, j = _lvira3d_optimize!(vg, ic, n_neigh, n0, work;
                                              igrid=igrid, tolir=tolir,
                                              niter=niter, wtype=wtype)

        if !(isfinite(nx) && isfinite(ny) && isfinite(nz) && isfinite(d) && isfinite(j))
            nx, ny, nz = n0
            d = _enforce_center_volume!(vg, ic, work, nx, ny, nz; igrid=igrid)
            if !isfinite(d)
                nx, ny, nz, d = _lsgir3d_cell!(vg, ic, work; pn=1.5, igrid=igrid)
            end
        end

        vg.xnormg[ic] = nx
        vg.ynormg[ic] = ny
        vg.znormg[ic] = nz
        vg.rholig[ic] = d
    end
    return vg
end
