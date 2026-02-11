# ---------------------------------------------------------------------------
#  User-defined test cases  (translation of usergvof.f)
#  – VARDEF, VELGRID, COMPDT, FUNC3DINIC / FUNC3DENDC / FUNCEXACTVOL
# ---------------------------------------------------------------------------

const π  = 3.141592653589793238462643383279
const π2 = 2π

# ═══════════════════════════════════════════════════════════════════════════
#  VARDEF – read parameters from a gvofvardef file
# ═══════════════════════════════════════════════════════════════════════════
"""
    vardef(filename::String="gvofvardef") -> VOFParams

Read a Fortran-style `gvofvardef` parameter file and return a `VOFParams`.
"""
function vardef(filename::String="gvofvardef")
    lines = readlines(filename)
    idx = 0
    _next() = (idx += 1; strip(lines[idx]))

    _next();     icase       = parse(Int, _next())
    _next();     irec        = parse(Int, _next())
    _next();     pn          = parse(Float64, _next())
    _next();     iw          = parse(Int, _next())
    _next();     parts       = split(_next())
                 niter       = parse(Int, parts[1])
                 tolir       = parse(Float64, parts[2])
    _next();     iadv        = parse(Int, _next())
    _next();     tcaden      = parse(Float64, _next())
    _next();     iprintplic  = parse(Int, _next()) == 1
    _next();     iprintiso   = parse(Int, _next()) == 1
    _next();     iprintvoxel = parse(Int, _next()) == 1
    _next();     nc          = parse(Int, _next())
    _next();     tolfr       = parse(Float64, _next())
    _next();     idt         = parse(Int, _next())
    _next();     dt          = parse(Float64, _next())
    _next();     cfl         = parse(Float64, _next())
    _next();     igrid       = parse(Int, _next())
    _next();     nx          = parse(Int, _next())
    _next();     ny          = parse(Int, _next())
    _next();     nz          = parse(Int, _next())

    xlast = ylast = zlast = 1.0
    if igrid == 1
        _next(); xlast = parse(Float64, _next())
        _next(); ylast = parse(Float64, _next())
        _next(); zlast = parse(Float64, _next())
    end

    timeend = if icase == 1;  4.0
         elseif icase == 2;  2π
         elseif icase == 3;  3.0
         elseif icase == 4;  8.0
         elseif icase == 5;  3.0
         else; 3.0; end

    return VOFParams(; icase, irec, pn, iw, niter, tolir, iadv, tcaden,
                       iprintplic, iprintiso, iprintvoxel, nc, tolfr,
                       idt, dt, cfl, igrid, nx, ny, nz,
                       xlast, ylast, zlast, timeend)
end

# ═══════════════════════════════════════════════════════════════════════════
#  VELGRID – evaluate test-case velocity fields
# ═══════════════════════════════════════════════════════════════════════════
"""
    velgrid!(vg::VOFGrid, icase::Int, timet0::Float64, timet::Float64)

Set face and node velocities for one of the five standard test cases.
Updates `vg.velface`, `vg.velnode` and `vg.velcmax`.
"""
function velgrid!(vg::VOFGrid, icase::Int, timet0::Float64, timet::Float64)
    g = vg.grid
    nface  = g.nface
    npoint = g.npoint
    tn = (timet0 + timet) / 2

    vg.velcmax .= 1e-16

    if icase == 1
        # Simple translation v = (0,0,1)
        for iface in 1:nface
            vg.velface[iface,1] = 0.0; vg.velface[iface,2] = 0.0; vg.velface[iface,3] = 1.0
        end
        for ip in 1:npoint
            vg.velnode[ip,1] = 0.0; vg.velnode[ip,2] = 0.0; vg.velnode[ip,3] = 1.0
        end
        vg.velcmax[3] = 1.0

    elseif icase == 2
        # XY rotation
        for iface in 1:nface
            x = vg.cface[iface,1] - 0.5; y = vg.cface[iface,2] - 0.5
            vg.velface[iface,1] = -y; vg.velface[iface,2] = x; vg.velface[iface,3] = 0.0
            vg.velcmax[1] = max(vg.velcmax[1], abs(vg.velface[iface,1]))
            vg.velcmax[2] = max(vg.velcmax[2], abs(vg.velface[iface,2]))
        end
        for ip in 1:npoint
            x = g.vnode[ip,1] - 0.5; y = g.vnode[ip,2] - 0.5
            vg.velnode[ip,1] = -y; vg.velnode[ip,2] = x; vg.velnode[ip,3] = 0.0
        end

    elseif icase == 3
        # 3D deformation (Enright et al. 2002)
        T = 3.0
        for iface in 1:nface
            xp, yp, zp = vg.cface[iface,1], vg.cface[iface,2], vg.cface[iface,3]
            vg.velface[iface,1] =  2*sin(π*xp)^2 * sin(π2*zp) * sin(π2*yp) * cos(π*tn/T)
            vg.velface[iface,2] = -sin(π*yp)^2 * sin(π2*xp) * sin(π2*zp) * cos(π*tn/T)
            vg.velface[iface,3] = -sin(π*zp)^2 * sin(π2*yp) * sin(π2*xp) * cos(π*tn/T)
            for i in 1:3; vg.velcmax[i] = max(vg.velcmax[i], abs(vg.velface[iface,i])); end
        end
        for ip in 1:npoint
            xp, yp, zp = g.vnode[ip,1], g.vnode[ip,2], g.vnode[ip,3]
            vg.velnode[ip,1] =  2*sin(π*xp)^2 * sin(π2*zp) * sin(π2*yp) * cos(π*tn/T)
            vg.velnode[ip,2] = -sin(π*yp)^2 * sin(π2*xp) * sin(π2*zp) * cos(π*tn/T)
            vg.velnode[ip,3] = -sin(π*zp)^2 * sin(π2*yp) * sin(π2*xp) * cos(π*tn/T)
        end

    elseif icase == 4
        # Vortex in a box XZ
        T = 8.0
        for iface in 1:nface
            xp, yp, zp = vg.cface[iface,1], vg.cface[iface,2], vg.cface[iface,3]
            vg.velface[iface,1] = -2*sin(π*xp)^2 * sin(π*zp) * cos(π*zp) * cos(π*tn/T)
            vg.velface[iface,2] =  0.0
            vg.velface[iface,3] =  2*sin(π*zp)^2 * sin(π*xp) * cos(π*xp) * cos(π*tn/T)
            vg.velcmax[1] = max(vg.velcmax[1], abs(vg.velface[iface,1]))
            vg.velcmax[3] = max(vg.velcmax[3], abs(vg.velface[iface,3]))
        end
        for ip in 1:npoint
            xp, yp, zp = g.vnode[ip,1], g.vnode[ip,2], g.vnode[ip,3]
            vg.velnode[ip,1] = -2*sin(π*xp)^2 * sin(π*zp) * cos(π*zp) * cos(π*tn/T)
            vg.velnode[ip,2] =  0.0
            vg.velnode[ip,3] =  2*sin(π*zp)^2 * sin(π*xp) * cos(π*xp) * cos(π*tn/T)
        end

    elseif icase == 5
        # 3D shearing (Liovic et al. 2006)
        T = 3.0
        for iface in 1:nface
            xp, yp, zp = vg.cface[iface,1], vg.cface[iface,2], vg.cface[iface,3]
            vg.velface[iface,1] = sin(π*xp)^2 * sin(π2*yp) * cos(π*tn/T)
            vg.velface[iface,2] = -sin(π*yp)^2 * sin(π2*xp) * cos(π*tn/T)
            vg.velface[iface,3] = (1.0 - 2*sqrt((xp-0.5)^2+(yp-0.5)^2))^2 * cos(π*tn/T)
            for i in 1:3; vg.velcmax[i] = max(vg.velcmax[i], abs(vg.velface[iface,i])); end
        end
        for ip in 1:npoint
            xp, yp, zp = g.vnode[ip,1], g.vnode[ip,2], g.vnode[ip,3]
            vg.velnode[ip,1] = sin(π*xp)^2 * sin(π2*yp) * cos(π*tn/T)
            vg.velnode[ip,2] = -sin(π*yp)^2 * sin(π2*xp) * cos(π*tn/T)
            vg.velnode[ip,3] = (1.0 - 2*sqrt((xp-0.5)^2+(yp-0.5)^2))^2 * cos(π*tn/T)
        end
    end

    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════
#  COMPDT – compute CFL-limited time step
# ═══════════════════════════════════════════════════════════════════════════
"""
    compdt!(vg::VOFGrid, params::VOFParams, timet::Float64) -> Float64

Compute the next time step `dt` from the CFL condition, with a 15 %
growth limiter and a final-time guard.  Returns the new dt and also
stores it in `params.dt` conceptually (caller must use the returned value).
"""
function compdt!(vg::VOFGrid, cfl::Float64, dt_prev::Float64,
                 timeend::Float64, timet::Float64)
    dt = cfl * min(vg.sizecmin[1] / vg.velcmax[1],
                   vg.sizecmin[2] / vg.velcmax[2],
                   vg.sizecmin[3] / vg.velcmax[3])
    dt_prev > 0 && (dt = min(dt, dt_prev * 1.15))
    timet + dt > timeend && (dt = timeend - timet)
    return dt
end

# ═══════════════════════════════════════════════════════════════════════════
#  Level-set shape functions for initial / end conditions
# ═══════════════════════════════════════════════════════════════════════════

# ── Initial shapes (φ > 0 inside) ──────────────────────────────────────
const _FUNC3DINIC = Dict{Int,Function}(
    1 => (x,y,z) -> -(  (x-0.5)^2 + (y-0.5)^2  + (z-0.5)^2  - 0.25^2 ),   # sphere r=0.25 at (0.5,0.5,0.5)
    2 => (x,y,z) -> -(  (x-0.5)^2 + (y-0.75)^2 + (z-0.5)^2  - 0.15^2 ),   # sphere r=0.15 at (0.5,0.75,0.5)
    3 => (x,y,z) -> -(  (x-0.35)^2+ (y-0.35)^2 + (z-0.35)^2 - 0.15^2 ),   # sphere r=0.15 at (0.35,0.35,0.35)
    4 => (x,y,z) -> -(  (x-0.5)^2 + (z-0.75)^2 - 0.15^2 ),                 # cylinder r=0.15 at (0.5,_,0.75)
    5 => (x,y,z) -> -(  (x-0.5)^2 + (y-0.75)^2 + (z-0.25)^2 - 0.15^2 ),   # sphere r=0.15 at (0.5,0.75,0.25)
    101 => (x,y,z) -> -( (x-0.525)^2 + (y-0.464)^2 + (z-0.516)^2 - 0.325^2 ), # sphere r=0.325
    102 => (x,y,z) -> begin   # hollow sphere
            φ1 = -( (x-0.525)^2+(y-0.464)^2+(z-0.516)^2 - 0.4^2 )
            φ2 =  ( (x-0.525)^2+(y-0.464)^2+(z-0.516)^2 - 0.2^2 )
            min(φ1, φ2)
        end,
    103 => (x,y,z) -> ( 0.1^2 - (0.2 - sqrt((x-0.525)^2+(y-0.464)^2))^2 - (z-0.516)^2 ),  # torus
)

# ── End-state shapes ───────────────────────────────────────────────────
const _FUNC3DENDC = Dict{Int,Function}(
    1 => (x,y,z) -> -( (x-0.5)^2 + (y-0.5)^2  + (z-4.5)^2  - 0.25^2 ),    # translated sphere
    2 => (x,y,z) -> -( (x-0.5)^2 + (y-0.75)^2 + (z-0.5)^2  - 0.15^2 ),
    3 => (x,y,z) -> -( (x-0.35)^2+ (y-0.35)^2 + (z-0.35)^2 - 0.15^2 ),
    4 => (x,y,z) -> -( (x-0.5)^2 + (z-0.75)^2 - 0.15^2 ),
    5 => (x,y,z) -> -( (x-0.5)^2 + (y-0.75)^2 + (z-0.25)^2 - 0.15^2 ),
)

# ── Exact volumes ──────────────────────────────────────────────────────
const _FUNCEXACTVOL = Dict{Int,Float64}(
    1 => 4π * 0.25^3 / 3,
    2 => 4π * 0.15^3 / 3,
    3 => 4π * 0.15^3 / 3,
    4 => π  * 0.15^2,              # per unit length (infinite cylinder)
    5 => 4π * 0.15^3 / 3,
    101 => 4π * 0.325^3 / 3,
    102 => 4π * (0.4^3 - 0.2^3) / 3,
    103 => 2 * π^2 * 0.1^2 * 0.2,
)

"""
    func3dinic(icase::Int) -> Function

Return the initial level-set function `(x,y,z) -> φ` for test case `icase`.
"""
func3dinic(icase::Int) = _FUNC3DINIC[icase]

"""
    func3dendc(icase::Int) -> Function

Return the end-state level-set function `(x,y,z) -> φ` for test case `icase`.
"""
func3dendc(icase::Int) = _FUNC3DENDC[icase]

"""
    funcexactvol(icase::Int) -> Float64

Return the analytical volume for test case `icase`.
"""
funcexactvol(icase::Int) = _FUNCEXACTVOL[icase]

# ═══════════════════════════════════════════════════════════════════════════
#  RECONSTRUCT! – dispatch to the appropriate reconstruction method
# ═══════════════════════════════════════════════════════════════════════════
"""
    reconstruct!(vg::VOFGrid, params::VOFParams)

Convenience function that calls the reconstruction method specified by
`params.irec` (1=CLCIR, 2=ELCIR, 3=LLCIR, 4=LSGIR, 5=SWIR, 6=LSFIR).
"""
function reconstruct!(vg::VOFGrid, params::VOFParams)
    irec  = params.irec
    igrid = params.igrid
    if irec == 1
        clcir!(vg; iw = params.iw, pn = params.pn, igrid = igrid)
    elseif irec == 2
        elcir!(vg; iw = params.iw, pn = params.pn, igrid = igrid)
    elseif irec == 3
        llcir!(vg; iw = params.iw, pn = params.pn, igrid = igrid)
    elseif irec == 4
        lsgir!(vg; pn = params.pn, igrid = igrid)
    elseif irec == 5
        swir!(vg; pn = params.pn, niter = params.niter,
              tolir = params.tolir, igrid = igrid)
    elseif irec == 6
        lsfir!(vg; pn = params.pn, niter = params.niter,
               tolir = params.tolir, igrid = igrid)
    else
        error("Unknown reconstruction method irec=$irec (expected 1–6)")
    end
    return nothing
end
