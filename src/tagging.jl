# ---------------------------------------------------------------------------
#  Tagging, initialisation and error computation
#  (TAGGRID, INITFGRID, RECERR from gvof.f)
# ---------------------------------------------------------------------------

# ═══════════════════════════════════════════════════════════════════════════
#  TAGGRID – IDW interpolation & three-level cell/face/node tagging
# ═══════════════════════════════════════════════════════════════════════════
"""
    taggrid!(vg::VOFGrid; tolfr=1e-12)

Classify every node, cell and face of the grid as *empty* (−1),
*interfacial* (0) or *full* (+1) and build the lists

* `vg.icint` – interfacial cells (those with 0 < f < 1)
* `vg.icadv` – cells that need advection
* `vg.isflu` – faces requiring a geometric flux computation

Nodal volume fractions `vg.frnod` are obtained by Inverse-Distance
Weighting (Shepard, p = 2) from the cell-centred values `vg.fractg`.
"""
function taggrid!(vg::VOFGrid; tolfr::Float64=1e-12)
    g = vg.grid
    ncell  = g.ncell
    nface  = g.nface
    npoint = g.npoint

    # ── 1) Node interpolation  (IDW, Shepard p = 2) ──────────────────
    for ip in 1:npoint
        vfnod = 0.0
        sw    = 0.0
        for ic in vg.icnode[ip]
            dist = sqrt((g.vnode[ip,1]-vg.ccell[ic,1])^2 +
                        (g.vnode[ip,2]-vg.ccell[ic,2])^2 +
                        (g.vnode[ip,3]-vg.ccell[ic,3])^2)
            w = 1.0 / dist^2
            vfnod += vg.fractg[ic] * w
            sw    += w
        end
        vg.frnod[ip] = vfnod / sw
        if vg.frnod[ip] < tolfr
            vg.iptag[ip] = -1
        elseif vg.frnod[ip] > (1.0 - tolfr)
            vg.iptag[ip] = 1
        else
            vg.iptag[ip] = 0
        end
    end

    # ── 2) Cell tagging ───────────────────────────────────────────────
    iadv = zeros(Int, ncell)
    for ic in 1:ncell
        if vg.fractg[ic] < tolfr
            vg.ictag[ic] = -1
        elseif vg.fractg[ic] > (1.0 - tolfr)
            vg.ictag[ic] = 1
        else
            vg.ictag[ic] = 0
        end
    end
    icint = Int[]
    for ic in 1:ncell
        vg.ictag[ic] == 0 && push!(icint, ic)
    end
    vg.icint = icint

    # ── 3) Face tagging ───────────────────────────────────────────────
    for iface in 1:nface
        ipf  = g.ipface[iface]
        nip  = length(ipf)
        isign = 1
        isum  = 0
        tagged = false
        for i in 1:nip
            ip = ipf[i]
            isign *= vg.iptag[ip]
            if isign == 0
                vg.istag[iface] = 0
                # Mark adjacent cells for advection
                ic = g.icface[iface,1]
                iadv[ic] == 0 && (iadv[ic] = 1)
                ic2 = g.icface[iface,2]
                ic2 != 0 && iadv[ic2] == 0 && (iadv[ic2] = 1)
                tagged = true
                break
            end
            isum += vg.iptag[ip]
        end
        if !tagged
            vg.istag[iface] = (isum == nip) ? 1 : -1
        end
    end

    isflu = Int[]
    for iface in 1:nface
        vg.istag[iface] == 0 && push!(isflu, iface)
    end
    vg.isflu = isflu

    icadv = Int[]
    for ic in 1:ncell
        iadv[ic] == 1 && push!(icadv, ic)
    end
    vg.icadv = icadv

    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════
#  INITFGRID – initialise volume fractions using a level-set function
# ═══════════════════════════════════════════════════════════════════════════
"""
    initfgrid!(vg::VOFGrid, func3d::Function; nc=10, tolfr=1e-12)

Initialise the volume fraction field `vg.fractg` by evaluating
`func3d(x, y, z)` (positive inside the material) over each cell using
`VOFTools.initf3d`.
"""
function initfgrid!(vg::VOFGrid, func3d::Function; nc::Int=10,
                    tolfr::Float64=1e-12)
    g   = vg.grid
    vft = 0.0
    for ic in 1:g.ncell
        poly = defcell(vg, ic)
        vf = initf3d(func3d, poly; nc=nc, tol=0.05)
        vf ≤ tolfr && (vf = 0.0)
        vf ≥ (1-tolfr) && (vf = 1.0)
        vg.fractg[ic] = vf
        vft += vf * vg.vcell[ic]
    end
    @info "Initialised volume: $vft"
    return vft
end

# ═══════════════════════════════════════════════════════════════════════════
#  RECERR – reconstruction error computation
# ═══════════════════════════════════════════════════════════════════════════
"""
    recerr(vg::VOFGrid, func3d::Function; nc=10) -> (erec, erinit, vexact)

Compute reconstruction and initialisation errors.

* `erec`   – symmetric-difference reconstruction error
* `erinit` – relative volume error |V_exact - V_init| / V_exact
* `vexact` – exact material volume (from `func3d`)
"""
function recerr(vg::VOFGrid, func3d::Function;
                nc::Int=10, vexact::Float64=NaN)
    g = vg.grid
    icint = vg.icint

    erec = 0.0
    for ind in eachindex(icint)
        ic = icint[ind]
        poly = defcell(vg, ic)
        xnc = vg.xnormg[ic]; ync = vg.ynormg[ic]; znc = vg.znormg[ic]
        c   = vg.rholig[ic]

        # Intersect PLIC with cell
        icontn, icontp = inte3d!(poly, c, xnc, ync, znc)
        vcut = toolv3d(poly)
        # Exact fraction of the PLIC-truncated cell
        vf = initf3d(func3d, poly; nc=nc, tol=0.05)
        erec += 2.0 * abs(vg.vcell[ic]*vg.fractg[ic] - vf*vcut)
    end

    # Total initialised volume
    vinit = 0.0
    for ic in 1:g.ncell
        vinit += vg.fractg[ic] * vg.vcell[ic]
    end

    if isnan(vexact)
        vex = vinit   # no reference → error = 0
    else
        vex = vexact
    end
    erinit = vex == 0 ? 0.0 : abs(vex - vinit) / vex

    return (erec=erec, erinit=erinit, vexact=vex)
end
