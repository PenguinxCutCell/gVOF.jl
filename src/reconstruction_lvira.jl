# ---------------------------------------------------------------------------
#  Canonical 2D ELVIRA / LVIRA reconstruction
# ---------------------------------------------------------------------------

const _ELVIRA_GRID_ERROR = "ELVIRA/LVIRA currently support only 2D uniform Cartesian grids with nz == 1."

@inline function _candidate_normal_from_slope_x(m::Float64)
    d = sqrt(1.0 + m * m)
    return (-m / d, 1.0 / d, 0.0)
end

@inline function _candidate_normal_from_slope_y(m::Float64)
    d = sqrt(1.0 + m * m)
    return (1.0 / d, -m / d, 0.0)
end

@inline function _normal_to_theta(nx::Float64, ny::Float64)
    return mod(atan(ny, nx), π)
end

@inline function _normal_from_theta(theta::Float64)
    return (cos(theta), sin(theta), 0.0)
end

@inline function _match_cartesian_index(values::Vector{Float64}, value::Float64;
                                        atol::Float64=1e-10,
                                        rtol::Float64=1e-10)
    idx = searchsortedfirst(values, value)
    best = 0
    besterr = Inf
    for j in (idx - 1, idx)
        (j < 1 || j > length(values)) && continue
        err = abs(values[j] - value)
        tol = max(atol, rtol * max(abs(values[j]), abs(value)))
        if err <= tol && err < besterr
            best = j
            besterr = err
        end
    end
    best == 0 && error(_ELVIRA_GRID_ERROR)
    return best
end

function _is_supported_elvira_grid(vg::VOFGrid, p=nothing)
    _ = p
    ncell = vg.grid.ncell
    ncell == 0 && return false

    dx = vg.boxcell[1,2] - vg.boxcell[1,1]
    dy = vg.boxcell[1,4] - vg.boxcell[1,3]
    dz = vg.boxcell[1,6] - vg.boxcell[1,5]
    if !(dx > 0.0 && dy > 0.0 && dz > 0.0)
        return false
    end

    for ic in 2:ncell
        if !isapprox(vg.boxcell[ic,2] - vg.boxcell[ic,1], dx; atol=1e-12, rtol=1e-10) ||
           !isapprox(vg.boxcell[ic,4] - vg.boxcell[ic,3], dy; atol=1e-12, rtol=1e-10) ||
           !isapprox(vg.boxcell[ic,6] - vg.boxcell[ic,5], dz; atol=1e-12, rtol=1e-10)
            return false
        end
    end

    z0 = vg.ccell[1,3]
    for ic in 2:ncell
        if !isapprox(vg.ccell[ic,3], z0; atol=1e-12, rtol=1e-10)
            return false
        end
    end

    xs = sort!(collect(unique(@view(vg.ccell[:,1]))))
    ys = sort!(collect(unique(@view(vg.ccell[:,2]))))
    zs = sort!(collect(unique(@view(vg.ccell[:,3]))))
    length(zs) == 1 || return false

    if length(xs) > 1
        dxs = diff(xs)
        for k in 2:length(dxs)
            isapprox(dxs[k], dxs[1]; atol=1e-12, rtol=1e-10) || return false
        end
    end
    if length(ys) > 1
        dys = diff(ys)
        for k in 2:length(dys)
            isapprox(dys[k], dys[1]; atol=1e-12, rtol=1e-10) || return false
        end
    end

    return true
end

function _build_cartesian_lookup(vg::VOFGrid, p=nothing)
    _ = p
    _is_supported_elvira_grid(vg, p) || error(_ELVIRA_GRID_ERROR)

    xs = sort!(collect(unique(@view(vg.ccell[:,1]))))
    ys = sort!(collect(unique(@view(vg.ccell[:,2]))))
    nx = length(xs)
    ny = length(ys)

    cell_of_ij = zeros(Int, nx, ny)
    ij_of_cell = Vector{NTuple{2,Int}}(undef, vg.grid.ncell)
    for ic in 1:vg.grid.ncell
        i = _match_cartesian_index(xs, vg.ccell[ic,1])
        j = _match_cartesian_index(ys, vg.ccell[ic,2])
        cell_of_ij[i, j] = ic
        ij_of_cell[ic] = (i, j)
    end

    any(==(0), cell_of_ij) && error(_ELVIRA_GRID_ERROR)
    return cell_of_ij, ij_of_cell
end

function _ensure_cartesian_lookup!(work::_ReconWork, vg::VOFGrid)
    cell_of_ij = work.cart_cell_of_ij
    ij_of_cell = work.cart_ij_of_cell
    need_build = size(cell_of_ij, 1) == 0 || size(cell_of_ij, 2) == 0 || length(ij_of_cell) != vg.grid.ncell
    if !need_build
        nx = size(cell_of_ij, 1)
        ny = size(cell_of_ij, 2)
        need_build = nx * ny != vg.grid.ncell
    end
    if need_build
        cell_of_ij, ij_of_cell = _build_cartesian_lookup(vg)
        work.cart_cell_of_ij = cell_of_ij
        work.cart_ij_of_cell = ij_of_cell
    end
    return work.cart_cell_of_ij, work.cart_ij_of_cell
end

function _full_3x3_stencil(cell_of_ij::Matrix{Int}, i::Int, j::Int)
    (i <= 1 || j <= 1 || i >= size(cell_of_ij, 1) || j >= size(cell_of_ij, 2)) && return nothing
    ids = (
        cell_of_ij[i - 1, j - 1], cell_of_ij[i, j - 1], cell_of_ij[i + 1, j - 1],
        cell_of_ij[i - 1, j],     cell_of_ij[i, j],     cell_of_ij[i + 1, j],
        cell_of_ij[i - 1, j + 1], cell_of_ij[i, j + 1], cell_of_ij[i + 1, j + 1],
    )
    any(==(0), ids) && return nothing
    return ids
end

function _fraction_for_cell_with_plane!(work::_ReconWork, vg::VOFGrid, ic::Int,
                                        nx::Float64, ny::Float64, nz::Float64,
                                        c::Float64)
    poly = defcell!(work.poly, work.ipg, vg, ic)
    inte3d!(poly, c, nx, ny, nz)
    vcut = toolv3d(poly)
    return vcut / vg.vcell[ic]
end

function _constrained_stencil_objective!(work::_ReconWork, vg::VOFGrid,
                                         stencil_ids::AbstractVector{Int},
                                         ic_center::Int,
                                         nx::Float64, ny::Float64, nz::Float64,
                                         c::Float64)
    _ = ic_center
    j = 0.0
    @inbounds for k in 1:9
        ic = stencil_ids[k]
        fhat = _fraction_for_cell_with_plane!(work, vg, ic, nx, ny, nz, c)
        df = fhat - vg.fractg[ic]
        j += work.stencil_weights[k] * df * df
    end
    return j
end

function _elvira_candidates_from_stencil(fracs3x3::AbstractVector{Float64})
    # gVOF's geometric clipping convention tracks the complementary phase
    # relative to canonical ELVIRA slope formulas.
    f11 = 1.0 - fracs3x3[1]; f21 = 1.0 - fracs3x3[2]; f31 = 1.0 - fracs3x3[3]
    f12 = 1.0 - fracs3x3[4]; f22 = 1.0 - fracs3x3[5]; f32 = 1.0 - fracs3x3[6]
    f13 = 1.0 - fracs3x3[7]; f23 = 1.0 - fracs3x3[8]; f33 = 1.0 - fracs3x3[9]

    c_left = f11 + f12 + f13
    c_mid  = f21 + f22 + f23
    c_right = f31 + f32 + f33

    r_bot = f11 + f21 + f31
    r_mid = f12 + f22 + f32
    r_top = f13 + f23 + f33

    mx_b = c_mid - c_left
    mx_c = 0.5 * (c_right - c_left)
    mx_f = c_right - c_mid

    my_b = r_mid - r_bot
    my_c = 0.5 * (r_top - r_bot)
    my_f = r_top - r_mid

    return (
        _candidate_normal_from_slope_x(mx_b),
        _candidate_normal_from_slope_x(mx_c),
        _candidate_normal_from_slope_x(mx_f),
        _candidate_normal_from_slope_y(my_b),
        _candidate_normal_from_slope_y(my_c),
        _candidate_normal_from_slope_y(my_f),
    )
end

function _local_lsgir2d!(work::_ReconWork, vg::VOFGrid, ic::Int; pn::Float64=1.5,
                         igrid::Int=1)
    poly = defcell!(work.poly, work.ipg, vg, ic)
    neigs = vg.ineigb[ic]
    fv = work.fv
    vx = work.vx
    vy = work.vy
    n = 0
    for ic2 in neigs
        n += 1
        fv[n] = vg.fractg[ic2]
        vx[n] = vg.ccell[ic2,1]
        vy[n] = vg.ccell[ic2,2]
    end
    n += 1
    fv[n] = vg.fractg[ic]
    vx[n] = vg.ccell[ic,1]
    vy[n] = vg.ccell[ic,2]

    xr = vx[n]
    yr = vy[n]
    fr = fv[n]
    ata11 = ata12 = ata22 = 0.0
    atb1 = atb2 = 0.0
    for i in 1:(n - 1)
        dx = vx[i] - xr
        dy = vy[i] - yr
        d = sqrt(dx * dx + dy * dy)
        d == 0.0 && continue
        w = 1.0 / d^pn
        a1 = w * dx
        a2 = w * dy
        b = w * (fv[i] - fr)
        ata11 += a1 * a1
        ata12 += a1 * a2
        ata22 += a2 * a2
        atb1 += a1 * b
        atb2 += a2 * b
    end

    det = ata11 * ata22 - ata12 * ata12
    if abs(det) ≤ 1e-30
        xn, yn, zn = 1.0, 0.0, 0.0
    else
        gx = (ata22 * atb1 - ata12 * atb2) / det
        gy = (ata11 * atb2 - ata12 * atb1) / det
        if !(isfinite(gx) && isfinite(gy))
            xn, yn, zn = 1.0, 0.0, 0.0
        else
            d = sqrt(gx * gx + gy * gy)
            if d > 0.0
                xn, yn, zn = gx / d, gy / d, 0.0
            else
                xn, yn, zn = 1.0, 0.0, 0.0
            end
        end
    end

    c = _enforce_volume(vg, ic, poly, xn, yn, zn; igrid=igrid)
    return xn, yn, zn, c
end

function _elvira_best_candidate!(work::_ReconWork, vg::VOFGrid,
                                 stencil_ids::AbstractVector{Int},
                                 fracs3x3::AbstractVector{Float64},
                                 ic_center::Int; igrid::Int=1)
    candidates = _elvira_candidates_from_stencil(fracs3x3)
    best_j = Inf
    best = nothing
    @inbounds for k in 1:6
        nx, ny, nz = candidates[k]
        if !(isfinite(nx) && isfinite(ny) && isfinite(nz))
            work.candidate_errors[k] = Inf
            continue
        end
        # Reload center-cell geometry for each candidate. Enforcement helpers
        # can mutate internal polyhedron buffers.
        poly = defcell!(work.poly, work.ipg, vg, ic_center)
        c = _enforce_volume(vg, ic_center, poly, nx, ny, nz; igrid=igrid)
        if !isfinite(c)
            work.candidate_errors[k] = Inf
            continue
        end
        j = _constrained_stencil_objective!(work, vg, stencil_ids, ic_center, nx, ny, nz, c)
        work.candidate_errors[k] = j
        if isfinite(j) && j < best_j
            best_j = j
            best = (nx, ny, nz, c, j)
        end
    end
    return best
end

function _objective_for_theta!(work::_ReconWork, vg::VOFGrid,
                               stencil_ids::AbstractVector{Int},
                               ic_center::Int, theta::Float64; igrid::Int=1)
    nx, ny, nz = _normal_from_theta(theta)
    poly = defcell!(work.poly, work.ipg, vg, ic_center)
    c = _enforce_volume(vg, ic_center, poly, nx, ny, nz; igrid=igrid)
    if !isfinite(c)
        return Inf, nx, ny, nz, c
    end
    j = _constrained_stencil_objective!(work, vg, stencil_ids, ic_center, nx, ny, nz, c)
    return j, nx, ny, nz, c
end

function _golden_section_minimize!(objective, a::Float64, b::Float64;
                                   tol::Float64=1e-10,
                                   niter::Int=64)
    ϕ = (sqrt(5.0) - 1.0) / 2.0
    c = b - ϕ * (b - a)
    d = a + ϕ * (b - a)
    fc = objective(c)
    fd = objective(d)
    for _ in 1:niter
        if abs(b - a) ≤ tol
            break
        end
        if fc ≤ fd
            b = d
            d = c
            fd = fc
            c = b - ϕ * (b - a)
            fc = objective(c)
        else
            a = c
            c = d
            fc = fd
            d = a + ϕ * (b - a)
            fd = objective(d)
        end
    end
    if fc ≤ fd
        return c, fc
    else
        return d, fd
    end
end

function _lvira_best_candidate!(work::_ReconWork, vg::VOFGrid,
                                stencil_ids::AbstractVector{Int},
                                fracs3x3::AbstractVector{Float64},
                                ic_center::Int; tolir::Float64=1e-10,
                                niter::Int=64, igrid::Int=1)
    elvira = _elvira_best_candidate!(work, vg, stencil_ids, fracs3x3, ic_center; igrid=igrid)
    elvira === nothing && return nothing

    nx0, ny0, nz0, c0, j0 = elvira
    theta0 = _normal_to_theta(nx0, ny0)
    local_width = π / 12.0

    objective(theta) = begin
        j, _, _, _, _ = _objective_for_theta!(work, vg, stencil_ids, ic_center, mod(theta, π); igrid=igrid)
        return j
    end

    local_j = Inf
    local_theta = theta0
    try
        local_theta, local_j = _golden_section_minimize!(objective, theta0 - local_width, theta0 + local_width;
                                                         tol=tolir, niter=niter)
    catch
        local_j = Inf
    end

    coarse_theta = theta0
    coarse_j = j0
    nsample = 32
    step = π / nsample
    for k in 0:(nsample - 1)
        theta = k * step
        j = objective(theta)
        if isfinite(j) && j < coarse_j
            coarse_j = j
            coarse_theta = theta
        end
    end

    best_theta = theta0
    best_j = j0
    if isfinite(local_j) && local_j < best_j
        best_theta = mod(local_theta, π)
        best_j = local_j
    end
    if isfinite(coarse_j) && coarse_j < best_j
        best_theta = coarse_theta
        best_j = coarse_j
    end

    if isfinite(coarse_j) && coarse_j < local_j
        try
            refined_theta, refined_j = _golden_section_minimize!(objective, coarse_theta - local_width,
                                                                coarse_theta + local_width;
                                                                tol=tolir, niter=niter)
            if isfinite(refined_j) && refined_j < best_j
                best_theta = mod(refined_theta, π)
                best_j = refined_j
            end
        catch
        end
    end

    j, nx, ny, nz, c = _objective_for_theta!(work, vg, stencil_ids, ic_center, best_theta; igrid=igrid)
    if !isfinite(j)
        return elvira
    end
    return (nx, ny, nz, c, j)
end

"""
    elvira!(vg::VOFGrid; igrid=1)

Canonical 2D ELVIRA reconstruction on uniform Cartesian grids with a single
z-layer (`nz == 1`). The interface position is enforced by exact center-cell
volume matching, and the same plane is evaluated across the full 3×3 stencil.

Cells on the boundary that do not have a complete 3×3 stencil fall back to a
local 2D least-squares gradient reconstruction for that cell only.
"""
function elvira!(vg::VOFGrid; igrid::Int=1)
    _is_supported_elvira_grid(vg) || error(_ELVIRA_GRID_ERROR)
    work = _get_recon_work(vg)
    cell_of_ij, ij_of_cell = _ensure_cartesian_lookup!(work, vg)

    stencil_ids = work.stencil_ids
    fracs3x3 = work.stencil_fracs

    for ic in vg.icint
        i, j = ij_of_cell[ic]
        stencil = _full_3x3_stencil(cell_of_ij, i, j)
        if stencil === nothing
            xn, yn, zn, c = _local_lsgir2d!(work, vg, ic; igrid=igrid)
            vg.xnormg[ic] = xn
            vg.ynormg[ic] = yn
            vg.znormg[ic] = zn
            vg.rholig[ic] = c
            continue
        end

        @inbounds for k in 1:9
            ic2 = stencil[k]
            stencil_ids[k] = ic2
            fracs3x3[k] = vg.fractg[ic2]
        end

        best = _elvira_best_candidate!(work, vg, stencil_ids, fracs3x3, ic; igrid=igrid)
        if best === nothing
            xn, yn, zn, c = _local_lsgir2d!(work, vg, ic; igrid=igrid)
        else
            xn, yn, zn, c, _ = best
        end

        vg.xnormg[ic] = xn
        vg.ynormg[ic] = yn
        vg.znormg[ic] = zn
        vg.rholig[ic] = c
    end
    return nothing
end

"""
    lvira!(vg::VOFGrid; igrid=1, tolir=1e-10, niter=64)

Canonical 2D LVIRA reconstruction on uniform Cartesian grids with a single
z-layer (`nz == 1`). A constrained 1D minimization is performed over the plane
angle after the center-cell equality constraint is enforced exactly.

If the nonlinear search fails, the implementation falls back to the ELVIRA
candidate for that cell. Boundary cells without a full 3×3 stencil fall back to
the local 2D least-squares reconstruction.
"""
function lvira!(vg::VOFGrid; igrid::Int=1, tolir::Float64=1e-10, niter::Int=64)
    _is_supported_elvira_grid(vg) || error(_ELVIRA_GRID_ERROR)
    work = _get_recon_work(vg)
    cell_of_ij, ij_of_cell = _ensure_cartesian_lookup!(work, vg)

    stencil_ids = work.stencil_ids
    fracs3x3 = work.stencil_fracs

    for ic in vg.icint
        i, j = ij_of_cell[ic]
        stencil = _full_3x3_stencil(cell_of_ij, i, j)
        if stencil === nothing
            xn, yn, zn, c = _local_lsgir2d!(work, vg, ic; igrid=igrid)
            vg.xnormg[ic] = xn
            vg.ynormg[ic] = yn
            vg.znormg[ic] = zn
            vg.rholig[ic] = c
            continue
        end

        @inbounds for k in 1:9
            ic2 = stencil[k]
            stencil_ids[k] = ic2
            fracs3x3[k] = vg.fractg[ic2]
        end

        best = _lvira_best_candidate!(work, vg, stencil_ids, fracs3x3, ic;
                                      tolir=tolir, niter=niter, igrid=igrid)
        if best === nothing
            xn, yn, zn, c = _local_lsgir2d!(work, vg, ic; igrid=igrid)
        else
            xn, yn, zn, c, _ = best
        end

        vg.xnormg[ic] = xn
        vg.ynormg[ic] = yn
        vg.znormg[ic] = zn
        vg.rholig[ic] = c
    end
    return nothing
end