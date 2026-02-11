# ---------------------------------------------------------------------------
#  Advection routines  (translation of EMFP, FMFP, NMFP, POLV,
#  FACEFLUX, CELLDFLUX, VOFADV from gvof.f)
# ---------------------------------------------------------------------------

# ═══════════════════════════════════════════════════════════════════════════
#  POLV – construct an adjustment polyhedron with a given volume
# ═══════════════════════════════════════════════════════════════════════════
"""
    polv(np::Int, vcor::Float64, base_verts::Matrix{Float64}) -> (ie, apex)

Construct an apex point that gives a pyramid of volume `vcor` above the
polygon defined by `base_verts[1:np,1:3]`.

Returns `(ie, apex)` where `ie=1` on success, `ie=0` if the polygon is
degenerate.
"""
function polv(np::Int, vcor::Float64, base_verts::Matrix{Float64})
    # Centroid of base polygon → vertex np+1
    cx = cy = cz = 0.0
    for ip in 1:np
        cx += base_verts[ip,1]; cy += base_verts[ip,2]; cz += base_verts[ip,3]
    end
    cx /= np; cy /= np; cz /= np

    # Newell-formula normal of the base polygon
    xn = yn = zn = 0.0
    for ip in 1:np
        ip2 = ip == np ? 1 : ip + 1
        xv1 = base_verts[ip2,1] - cx;  yv1 = base_verts[ip2,2] - cy;  zv1 = base_verts[ip2,3] - cz
        xv2 = base_verts[ip,1]  - cx;  yv2 = base_verts[ip,2]  - cy;  zv2 = base_verts[ip,3]  - cz
        xn += yv1*zv2 - zv1*yv2
        yn += zv1*xv2 - xv1*zv2
        zn += xv1*yv2 - yv1*xv2
    end
    if xn == 0.0 && yn == 0.0 && zn == 0.0
        return (0, (cx, cy, cz))
    end
    xnc = xn / np; ync = yn / np; znc = zn / np
    c = 6.0 * vcor / (np * (xnc^2 + ync^2 + znc^2))
    apex = (cx + c*xnc, cy + c*ync, cz + c*znc)
    return (1, apex)
end

# ═══════════════════════════════════════════════════════════════════════════
#  Helper – compute face normal vectors for a raw polyhedron
# ═══════════════════════════════════════════════════════════════════════════
function _compute_normals!(xns, yns, zns, ipv, nipv, nts, vertp)
    for is in 1:nts
        ip1 = ipv[is,1]; ip2 = ipv[is,2]; ip3 = ipv[is,3]
        xv1 = vertp[ip2,1]-vertp[ip1,1]; yv1 = vertp[ip2,2]-vertp[ip1,2]; zv1 = vertp[ip2,3]-vertp[ip1,3]
        xv2 = vertp[ip3,1]-vertp[ip2,1]; yv2 = vertp[ip3,2]-vertp[ip2,2]; zv2 = vertp[ip3,3]-vertp[ip2,3]
        xm = yv1*zv2-zv1*yv2; ym = zv1*xv2-xv1*zv2; zm = xv1*yv2-yv1*xv2
        amod = sqrt(xm^2+ym^2+zm^2)
        if amod == 0
            xns[is] = yns[is] = zns[is] = 0.0
        else
            xns[is] = xm/amod; yns[is] = ym/amod; zns[is] = zm/amod
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════
#  EMFP – Edge-Matched Flux Polyhedron
# ═══════════════════════════════════════════════════════════════════════════
"""
    emfp!(dt, nivface, vertp, vel, v) -> (ipv, nipv, nts, ntp, vertp_out, xns, yns, zns, vfactor)

Construct the edge-matched flux polyhedron.

* `vertp[1:nivface,1:3]` – face vertex coordinates
* `vel[1:nivface,1:3]`   – per-vertex velocity
* `v`                     – signed face-flux volume

Returns arrays defining the polyhedron and a volume-correction factor.
"""
function emfp!(dt::Float64, nivface::Int, vertp_in::Matrix{Float64},
               vel::Matrix{Float64}, v::Float64)
    nf = nivface
    nts = nf*(4+1) + 1   # faces
    ntp = 3*nf + 1       # vertices
    ntv = ntp

    # Allocate with VOFTools default sizes so inte3d!/newpol3d! have room
    ns_max = VOFTools.NS_DEFAULT
    nv_max = VOFTools.NV_DEFAULT
    vertp = zeros(nv_max, 3)
    ipv   = zeros(Int, ns_max, nv_max)
    nipv  = zeros(Int, ns_max)
    xns   = zeros(ns_max); yns = zeros(ns_max); zns = zeros(ns_max)

    # Copy front face
    for ip in 1:nf, i in 1:3
        vertp[ip,i] = vertp_in[ip,i]
    end

    # Streak-lines → back vertices
    for ip in 1:nf, i in 1:3
        vertp[nf+ip,i] = vertp[ip,i] - vel[ip,i]*dt
    end

    # Geometric centres of each lateral quad (midpoints)
    vertp[3*nf+1,:] .= 0.0
    for ip in 1:nf
        ip1 = ip == nf ? 1 : ip + 1
        for i in 1:3
            vertp[2*nf+ip,i] = (vertp[ip,i]+vertp[nf+ip,i]+vertp[ip1,i]+vertp[nf+ip1,i])/4.0
            vertp[3*nf+1,i] += vertp[nf+ip,i]/nf
        end
    end

    # Front face
    nipv[1] = nf
    for ip in 1:nf; ipv[1,ip] = ip; end

    # Lateral faces (4 triangles per edge)
    for ip in 1:nf
        ip1 = ip == nf ? 1 : ip + 1
        is = (ip-1)*4+2
        nipv[is] = 3; ipv[is,1] = ip; ipv[is,2] = nf+ip; ipv[is,3] = 2*nf+ip

        is = (ip-1)*4+3
        nipv[is] = 3; ipv[is,1] = ip1; ipv[is,2] = ip; ipv[is,3] = 2*nf+ip

        is = (ip-1)*4+4
        nipv[is] = 3; ipv[is,1] = nf+ip1; ipv[is,2] = ip1; ipv[is,3] = 2*nf+ip

        is = (ip-1)*4+5
        nipv[is] = 3; ipv[is,1] = nf+ip; ipv[is,2] = nf+ip1; ipv[is,3] = 2*nf+ip
    end

    # Back face (triangles radiating from geometric centre)
    for ip in 1:nf
        ip1 = ip == nf ? 1 : ip + 1
        is = 1+4*nf+ip
        nipv[is] = 3; ipv[is,1] = nf+ip1; ipv[is,2] = nf+ip; ipv[is,3] = 3*nf+1
    end

    # Face normal vectors
    _compute_normals!(xns, yns, zns, ipv, nipv, nts, vertp)

    # Volume adjustment
    vini = toolv3d(ipv, nipv, nts, vertp, xns, yns, zns)
    vcor = v - vini
    base = zeros(nf+2, 3)
    for ip in 1:nf, i in 1:3
        base[ip,i] = vertp[nf+ip,i]
    end
    ie, apex = polv(nf, vcor, base)

    vfactor = 1.0
    if ie == 0
        vfactor = vini == 0 ? 1.0 : v / vini
    else
        vertp[ntv,1] = apex[1]; vertp[ntv,2] = apex[2]; vertp[ntv,3] = apex[3]
        # Re-compute normals for back faces
        for is in (nts-nf+1):nts
            _compute_normals!(xns, yns, zns, ipv, nipv, is:is, vertp)
        end
    end

    return (ipv, nipv, nts, ntp, vertp, xns, yns, zns, vfactor)
end

# Single-face normal update helper
function _compute_normals!(xns, yns, zns, ipv, nipv, rng::UnitRange{Int}, vertp)
    for is in rng
        ip1 = ipv[is,1]; ip2 = ipv[is,2]; ip3 = ipv[is,3]
        xv1 = vertp[ip2,1]-vertp[ip1,1]; yv1 = vertp[ip2,2]-vertp[ip1,2]; zv1 = vertp[ip2,3]-vertp[ip1,3]
        xv2 = vertp[ip3,1]-vertp[ip2,1]; yv2 = vertp[ip3,2]-vertp[ip2,2]; zv2 = vertp[ip3,3]-vertp[ip2,3]
        xm = yv1*zv2-zv1*yv2; ym = zv1*xv2-xv1*zv2; zm = xv1*yv2-yv1*xv2
        amod = sqrt(xm^2+ym^2+zm^2)
        if amod == 0
            xns[is] = yns[is] = zns[is] = 0.0
        else
            xns[is] = xm/amod; yns[is] = ym/amod; zns[is] = zm/amod
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════
#  FMFP – Face-Matched Flux Polyhedron
# ═══════════════════════════════════════════════════════════════════════════
"""
    fmfp!(dt, nivface, vertp_in, vel, v) -> (ipv, nipv, nts, ntp, vertp, xns, yns, zns, vfactor)

Construct the face-matched flux polyhedron.
"""
function fmfp!(dt::Float64, nivface::Int, vertp_in::Matrix{Float64},
               vel::Matrix{Float64}, v::Float64)
    nf = nivface
    nts = nf*2 + 1
    ntp = 2*nf + 1
    ntv = ntp

    # Allocate with VOFTools default sizes so inte3d!/newpol3d! have room
    ns_max = VOFTools.NS_DEFAULT
    nv_max = VOFTools.NV_DEFAULT
    vertp = zeros(nv_max, 3)
    ipv   = zeros(Int, ns_max, nv_max)
    nipv  = zeros(Int, ns_max)
    xns   = zeros(ns_max); yns = zeros(ns_max); zns = zeros(ns_max)

    for ip in 1:nf, i in 1:3
        vertp[ip,i] = vertp_in[ip,i]
    end

    # Orientation of lateral faces from edge × velocity
    edge = zeros(nf, 3)
    vele = zeros(nf, 3)
    for ip in 1:nf
        ip2 = ip == nf ? 1 : ip + 1
        for i in 1:3; vele[ip,i] = (vel[ip,i]+vel[ip2,i])/2; end
        xv = vertp[ip2,1]-vertp[ip,1]; yv = vertp[ip2,2]-vertp[ip,2]; zv = vertp[ip2,3]-vertp[ip,3]
        xn = yv*vele[ip,3]-zv*vele[ip,2]
        yn = zv*vele[ip,1]-xv*vele[ip,3]
        zn = xv*vele[ip,2]-yv*vele[ip,1]
        dmod = sqrt(xn^2+yn^2+zn^2)
        if dmod != 0
            xns[ip+1] = xn/dmod; yns[ip+1] = yn/dmod; zns[ip+1] = zn/dmod
        end
    end

    # Edge directions + new vertex positions
    for ip in 1:nf
        is0 = ip == 1 ? nf+1 : ip
        is = ip + 1
        edge[ip,1] = yns[is0]*zns[is]-zns[is0]*yns[is]
        edge[ip,2] = zns[is0]*xns[is]-xns[is0]*zns[is]
        edge[ip,3] = xns[is0]*yns[is]-yns[is0]*xns[is]
        dmod = sqrt(edge[ip,1]^2+edge[ip,2]^2+edge[ip,3]^2)
        if dmod != 0
            edge[ip,:] ./= dmod
        end
        velmod = dt*(vel[ip,1]*edge[ip,1]+vel[ip,2]*edge[ip,2]+vel[ip,3]*edge[ip,3])
        for i in 1:3
            vertp[nf+ip,i] = vertp[ip,i] - velmod*edge[ip,i]
        end
    end

    # Back face centre
    vertp[2*nf+1,:] .= 0.0
    for ip in 1:nf, i in 1:3
        vertp[2*nf+1,i] += vertp[nf+ip,i]/nf
    end

    # Front face
    nipv[1] = nf
    for ip in 1:nf; ipv[1,ip] = ip; end

    # Lateral quad faces
    for ip in 1:nf
        ip1 = ip == nf ? 1 : ip + 1
        is = ip + 1
        nipv[is] = 4
        ipv[is,1] = ip; ipv[is,2] = nf+ip; ipv[is,3] = nf+ip1; ipv[is,4] = ip1
    end

    # Back face triangles
    for ip in 1:nf
        ip1 = ip == nf ? 1 : ip + 1
        is = 1+nf+ip
        nipv[is] = 3; ipv[is,1] = nf+ip1; ipv[is,2] = nf+ip; ipv[is,3] = 2*nf+1
    end

    # Normals for back faces
    for is in (nf+2):nts
        _compute_normals!(xns, yns, zns, ipv, nipv, is:is, vertp)
    end

    # Volume adjustment
    vini = toolv3d(ipv, nipv, nts, vertp, xns, yns, zns)
    vcor = v - vini
    base = zeros(nf+2, 3)
    for ip in 1:nf, i in 1:3
        base[ip,i] = vertp[nf+ip,i]
    end
    ie, apex = polv(nf, vcor, base)

    vfactor = 1.0
    if ie == 0
        vfactor = vini == 0 ? 1.0 : v / vini
    else
        vertp[ntv,1] = apex[1]; vertp[ntv,2] = apex[2]; vertp[ntv,3] = apex[3]
        for is in (nts-nf+1):nts
            _compute_normals!(xns, yns, zns, ipv, nipv, is:is, vertp)
        end
    end

    return (ipv, nipv, nts, ntp, vertp, xns, yns, zns, vfactor)
end

# ═══════════════════════════════════════════════════════════════════════════
#  NMFP – Non-Matched Flux Polyhedron (Rider–Kothe extended to 3-D)
# ═══════════════════════════════════════════════════════════════════════════
"""
    nmfp!(dt, nivface, vertp_in, velface, v) -> (ipv, nipv, nts, ntp, vertp, xns, yns, zns, vfactor)

Construct the non-matched flux polyhedron (uniform back-translation).
`velface` is a 3-vector (face-centred velocity).
"""
function nmfp!(dt::Float64, nivface::Int, vertp_in::Matrix{Float64},
               velface::AbstractVector{Float64}, v::Float64)
    nf = nivface
    nts = nf + 2
    ntp = 2*nf
    ntv = ntp

    # Allocate with VOFTools default sizes so inte3d!/newpol3d! have room
    ns_max = VOFTools.NS_DEFAULT
    nv_max = VOFTools.NV_DEFAULT
    vertp = zeros(nv_max, 3)
    ipv   = zeros(Int, ns_max, nv_max)
    nipv  = zeros(Int, ns_max)
    xns   = zeros(ns_max); yns = zeros(ns_max); zns = zeros(ns_max)

    for ip in 1:nf, i in 1:3
        vertp[ip,i] = vertp_in[ip,i]
    end

    # Streak-lines (uniform velocity)
    for ip in 1:nf, i in 1:3
        vertp[nf+ip,i] = vertp[ip,i] - velface[i]*dt
    end

    # Front face
    nipv[1] = nf
    for ip in 1:nf; ipv[1,ip] = ip; end

    # Lateral quad faces
    for ip in 1:nf
        ip1 = ip == nf ? 1 : ip + 1
        is = ip + 1
        nipv[is] = 4
        ipv[is,1] = ip; ipv[is,2] = nf+ip; ipv[is,3] = nf+ip1; ipv[is,4] = ip1
    end

    # Back face
    is = nts; nipv[is] = nf
    for ip in 1:nf; ipv[is,ip] = ntp-ip+1; end

    # Face normal vectors
    _compute_normals!(xns, yns, zns, ipv, nipv, nts, vertp)

    # Volume factor
    vini = toolv3d(ipv, nipv, nts, vertp, xns, yns, zns)
    vfactor = vini == 0.0 ? 1.0 : v / vini

    return (ipv, nipv, nts, ntp, vertp, xns, yns, zns, vfactor)
end

# ═══════════════════════════════════════════════════════════════════════════
#  CELLDFLUX – cell decomposition for flux computation
# ═══════════════════════════════════════════════════════════════════════════
"""
    celldflux(vg, ic, igrid, xmin, xmax, ymin, ymax, zmin, zmax)
        -> Vector of (isign, faces)

For `igrid==3` (non-convex cells), decomposes the cell into convex sub-cells.
Otherwise returns a single group containing all cell faces.

Each element is `(isign, xnf, ynf, znf, rhof)` vectors.
"""
function celldflux(vg::VOFGrid, ic::Int, igrid::Int,
                   xmin::Float64, xmax::Float64,
                   ymin::Float64, ymax::Float64,
                   zmin::Float64, zmax::Float64)
    g = vg.grid
    if igrid == 3
        # Decompose non-convex cell into sub-tetrahedra (face + centroid)
        result = NamedTuple{(:isign,:xnf,:ynf,:znf,:rhof),
                            Tuple{Int,Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64}}}[]
        for is in 1:length(g.iscell[ic])
            iface = g.iscell[ic][is]
            nip = length(g.ipface[iface])

            # Bounding-box check
            xmin2 = ymin2 = zmin2 =  1e16
            xmax2 = ymax2 = zmax2 = -1e16
            xmin2 = min(xmin2, vg.ccell[ic,1]); xmax2 = max(xmax2, vg.ccell[ic,1])
            ymin2 = min(ymin2, vg.ccell[ic,2]); ymax2 = max(ymax2, vg.ccell[ic,2])
            zmin2 = min(zmin2, vg.ccell[ic,3]); zmax2 = max(zmax2, vg.ccell[ic,3])
            for idx in 1:nip
                ip1 = g.ipface[iface][idx]
                xmin2 = min(xmin2, g.vnode[ip1,1]); xmax2 = max(xmax2, g.vnode[ip1,1])
                ymin2 = min(ymin2, g.vnode[ip1,2]); ymax2 = max(ymax2, g.vnode[ip1,2])
                zmin2 = min(zmin2, g.vnode[ip1,3]); zmax2 = max(zmax2, g.vnode[ip1,3])
            end
            (xmax2 ≤ xmin || xmin2 ≥ xmax ||
             ymax2 ≤ ymin || ymin2 ≥ ymax ||
             zmax2 ≤ zmin || zmin2 ≥ zmax) && continue

            # Build half-plane decomposition
            niscelld = nip + 1
            xnf = zeros(niscelld); ynf = zeros(niscelld); znf = zeros(niscelld)
            rhof = zeros(niscelld)

            # First plane: the face itself (oriented towards cell centroid)
            if g.icface[iface,1] == ic
                xnf[1] = -vg.xnface[iface]; ynf[1] = -vg.ynface[iface]; znf[1] = -vg.znface[iface]
            else
                xnf[1] = vg.xnface[iface]; ynf[1] = vg.ynface[iface]; znf[1] = vg.znface[iface]
            end
            ip0 = g.ipface[iface][1]
            rhof[1] = -(xnf[1]*g.vnode[ip0,1]+ynf[1]*g.vnode[ip0,2]+znf[1]*g.vnode[ip0,3])

            # Remaining planes: edge + centroid
            if g.icface[iface,1] == ic
                order = 1:nip; inc = 1
            else
                order = nip:-1:1; inc = -1
            end
            is2 = 1
            for idx in order
                is2 += 1
                ip1 = g.ipface[iface][idx]
                idx2 = idx + inc; idx2 == 0 && (idx2 = nip); idx2 == nip+1 && (idx2 = 1)
                ip2 = g.ipface[iface][idx2]
                xv1 = g.vnode[ip2,1]-g.vnode[ip1,1]; yv1 = g.vnode[ip2,2]-g.vnode[ip1,2]; zv1 = g.vnode[ip2,3]-g.vnode[ip1,3]
                xv2 = vg.ccell[ic,1]-g.vnode[ip2,1]; yv2 = vg.ccell[ic,2]-g.vnode[ip2,2]; zv2 = vg.ccell[ic,3]-g.vnode[ip2,3]
                xn = yv1*zv2-zv1*yv2; yn = zv1*xv2-xv1*zv2; zn = xv1*yv2-yv1*xv2
                dmod = sqrt(xn^2+yn^2+zn^2)
                if dmod != 0
                    xnf[is2] = xn/dmod; ynf[is2] = yn/dmod; znf[is2] = zn/dmod
                end
                rhof[is2] = -(xnf[is2]*g.vnode[ip1,1]+ynf[is2]*g.vnode[ip1,2]+znf[is2]*g.vnode[ip1,3])
            end

            # Sign check: cell centroid should be on positive side of first plane
            isign = 1
            if rhof[1]+xnf[1]*vg.ccell[ic,1]+ynf[1]*vg.ccell[ic,2]+znf[1]*vg.ccell[ic,3] < 0
                isign = -1
                xnf .= -xnf; ynf .= -ynf; znf .= -znf; rhof .= -rhof
            end

            push!(result, (isign=isign, xnf=xnf, ynf=ynf, znf=znf, rhof=rhof))
        end
        return result
    else
        # Convex cell: single group of face half-planes
        nsc = length(g.iscell[ic])
        xnf = zeros(nsc); ynf = zeros(nsc); znf = zeros(nsc); rhof = zeros(nsc)
        for is in 1:nsc
            iface = g.iscell[ic][is]
            if g.icface[iface,1] == ic
                xnf[is] = -vg.xnface[iface]; ynf[is] = -vg.ynface[iface]; znf[is] = -vg.znface[iface]
            else
                xnf[is] = vg.xnface[iface]; ynf[is] = vg.ynface[iface]; znf[is] = vg.znface[iface]
            end
            ip0 = g.ipface[iface][1]
            rhof[is] = -(xnf[is]*g.vnode[ip0,1]+ynf[is]*g.vnode[ip0,2]+znf[is]*g.vnode[ip0,3])
        end
        return [(isign=1, xnf=xnf, ynf=ynf, znf=znf, rhof=rhof)]
    end
end

# ═══════════════════════════════════════════════════════════════════════════
#  FACEFLUX – face flux computation
# ═══════════════════════════════════════════════════════════════════════════
"""
    faceflux!(vg::VOFGrid; dt, iadv=1, igrid=1, tolfr=1e-12)

Compute the face volumetric fluxes `vg.fface` for every face in
`vg.isflu`.  Uses the selected advection method (EMFP/FMFP/NMFP).
"""
function faceflux!(vg::VOFGrid; dt::Float64, iadv::Int=1,
                   igrid::Int=1, tolfr::Float64=1e-12)
    g     = vg.grid
    nface = g.nface
    isflu = vg.isflu

    # Pre-compute signed volumes for all faces
    for iface in 1:nface
        vg.volpol[iface] = vg.aface[iface] * dt * (
            vg.velface[iface,1]*vg.xnface[iface] +
            vg.velface[iface,2]*vg.ynface[iface] +
            vg.velface[iface,3]*vg.znface[iface])
        if vg.istag[iface] == -1
            vg.fface[iface] = 0.0
        elseif vg.istag[iface] == 1
            vg.fface[iface] = vg.volpol[iface]
        end
    end

    for ind in eachindex(isflu)
        iface = isflu[ind]
        nivface = length(g.ipface[iface])

        # Gather face vertex coords and velocities
        face_verts = zeros(nivface, 3)
        face_vel   = zeros(nivface, 3)
        for i in 1:nivface
            ip = g.ipface[iface][i]
            face_verts[i,:] .= g.vnode[ip,:]
            face_vel[i,:]   .= vg.velnode[ip,:]
        end

        v = vg.volpol[iface]

        # Build flux polyhedron
        if iadv == 1
            ipv, nipv, nts, ntp, vertp, xns, yns, zns, vfactor =
                emfp!(dt, nivface, face_verts, face_vel, v)
        elseif iadv == 2
            ipv, nipv, nts, ntp, vertp, xns, yns, zns, vfactor =
                fmfp!(dt, nivface, face_verts, face_vel, v)
        else
            ipv, nipv, nts, ntp, vertp, xns, yns, zns, vfactor =
                nmfp!(dt, nivface, face_verts, @view(vg.velface[iface,:]), v)
        end

        # Bounding box of flux polyhedron
        fp_xmin = fp_ymin = fp_zmin =  1e16
        fp_xmax = fp_ymax = fp_zmax = -1e16
        for ip in 1:ntp
            fp_xmin = min(fp_xmin, vertp[ip,1]); fp_xmax = max(fp_xmax, vertp[ip,1])
            fp_ymin = min(fp_ymin, vertp[ip,2]); fp_ymax = max(fp_ymax, vertp[ip,2])
            fp_zmin = min(fp_zmin, vertp[ip,3]); fp_zmax = max(fp_zmax, vertp[ip,3])
        end

        # Reference cell (upstream of the face)
        dsign = vg.xnface[iface]*vg.velface[iface,1] +
                vg.ynface[iface]*vg.velface[iface,2] +
                vg.znface[iface]*vg.velface[iface,3]
        if dsign ≥ 0 || g.icface[iface,2] == 0
            icref = g.icface[iface,1]
        else
            icref = g.icface[iface,2]
        end

        # Build candidate cell list via bounding-box checks
        icflux0 = Int[]   # cells with f ≈ 0
        icflux1 = Int[]   # cells with f ≈ 1
        icfluxi = Int[]   # interface cells
        tot_xmin = tot_ymin = tot_zmin =  1e16
        tot_xmax = tot_ymax = tot_zmax = -1e16

        candidates = push!(copy(vg.ineigb[icref]), icref)
        for ic in candidates
            cx1 = vg.boxcell[ic,1]; cx2 = vg.boxcell[ic,2]
            cy1 = vg.boxcell[ic,3]; cy2 = vg.boxcell[ic,4]
            cz1 = vg.boxcell[ic,5]; cz2 = vg.boxcell[ic,6]
            tot_xmin = min(tot_xmin, cx1); tot_xmax = max(tot_xmax, cx2)
            tot_ymin = min(tot_ymin, cy1); tot_ymax = max(tot_ymax, cy2)
            tot_zmin = min(tot_zmin, cz1); tot_zmax = max(tot_zmax, cz2)
            (cx2 ≤ fp_xmin || cx1 ≥ fp_xmax ||
             cy2 ≤ fp_ymin || cy1 ≥ fp_ymax ||
             cz2 ≤ fp_zmin || cz1 ≥ fp_zmax) && continue
            if vg.fractg[ic] < tolfr
                push!(icflux0, ic)
            elseif vg.fractg[ic] > (1-tolfr)
                push!(icflux1, ic)
            else
                push!(icfluxi, ic)
                push!(icflux0, ic)
                push!(icflux1, ic)
            end
        end

        # Compute volumetric flux
        vflu = 0.0
        n0 = length(icflux0) - length(icfluxi)
        n1 = length(icflux1) - length(icfluxi)
        ni = length(icfluxi)

        if (ni + n0 + n1) != 0
            if ni == 0 && n0 == 0
                vflu = v
            else
                # Choose complementary direction that minimises work
                if n0 < n1 && fp_xmin ≥ tot_xmin && fp_xmax ≤ tot_xmax &&
                   fp_ymin ≥ tot_ymin && fp_ymax ≤ tot_ymax &&
                   fp_zmin ≥ tot_zmin && fp_zmax ≤ tot_zmax
                    signflux = -1.0
                    icflux = icflux0
                else
                    signflux = 1.0
                    icflux = icflux1
                end

                for ic in icflux
                    # Copy current flux polyhedron for intersection
                    ipv0  = copy(ipv); nipv0 = copy(nipv)
                    nts0  = nts; ntp0 = ntp
                    vertp0 = copy(vertp); xns0 = copy(xns); yns0 = copy(yns); zns0 = copy(zns)

                    # Intersect with PLIC interface
                    if vg.fractg[ic] > tolfr && vg.fractg[ic] < (1-tolfr)
                        xnc = vg.xnormg[ic]*signflux
                        ync = vg.ynormg[ic]*signflux
                        znc = vg.znormg[ic]*signflux
                        c_plic = vg.rholig[ic]*signflux

                        poly0 = Polyhedron3D(vertp0, ipv0, nipv0, xns0, yns0, zns0,
                                            nts0, ntp0, ntp0)
                        icontn, icontp = inte3d!(poly0, c_plic, xnc, ync, znc)
                        icontp == 0 && @goto next_cell
                        nts0 = poly0.nts; ntp0 = poly0.ntp
                        vertp0 = poly0.vertp; ipv0 = poly0.ipv; nipv0 = poly0.nipv
                        xns0 = poly0.xns; yns0 = poly0.yns; zns0 = poly0.zns
                    end

                    # Intersect with cell faces
                    decomp = celldflux(vg, ic, igrid, fp_xmin, fp_xmax,
                                       fp_ymin, fp_ymax, fp_zmin, fp_zmax)

                    ipv_save = nothing
                    if igrid == 3
                        ipv_save = (copy(ipv0), copy(nipv0), nts0, ntp0,
                                    copy(vertp0), copy(xns0), copy(yns0), copy(zns0))
                    end

                    for (di, d) in enumerate(decomp)
                        if igrid == 3 && di > 1
                            # Restore from save
                            ipv0, nipv0, nts0, ntp0, vertp0, xns0, yns0, zns0 = ipv_save
                            ipv0 = copy(ipv0); nipv0 = copy(nipv0)
                            vertp0 = copy(vertp0); xns0 = copy(xns0); yns0 = copy(yns0); zns0 = copy(zns0)
                        end

                        cont = false
                        for i in eachindex(d.xnf)
                            poly0 = Polyhedron3D(vertp0, ipv0, nipv0, xns0, yns0, zns0,
                                                nts0, ntp0, ntp0)
                            icontn, icontp = inte3d!(poly0, d.rhof[i],
                                                     d.xnf[i], d.ynf[i], d.znf[i])
                            nts0 = poly0.nts; ntp0 = poly0.ntp
                            vertp0 = poly0.vertp; ipv0 = poly0.ipv; nipv0 = poly0.nipv
                            xns0 = poly0.xns; yns0 = poly0.yns; zns0 = poly0.zns
                            if icontp == 0
                                cont = true; break
                            end
                        end
                        cont && continue

                        # Compute truncated volume
                        vol_trunc = toolv3d(ipv0, nipv0, nts0, vertp0, xns0, yns0, zns0)
                        vflu += vol_trunc * d.isign
                    end
                    @label next_cell
                end
            end
        end

        vg.fface[iface] = vflu * vfactor
        if ni > 0 && signflux < 0
            vg.fface[iface] = (vg.volpol[iface] - vflu) * vfactor
        end
    end
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════
#  VOFADV – VOF advection with implicit divergence correction
# ═══════════════════════════════════════════════════════════════════════════
"""
    vofadv!(vg::VOFGrid; tolfr=1e-12)

Advance the volume fraction field using the face fluxes `vg.fface`.
Uses implicit divergence correction.

Updates `vg.fractg` and accumulates boundedness errors in `vg.ebound`.
"""
function vofadv!(vg::VOFGrid; tolfr::Float64=1e-12)
    g     = vg.grid
    ncell = g.ncell
    icadv = vg.icadv

    fract0 = copy(vg.fractg)
    funder = 0.0
    fover  = 0.0

    for ind in eachindex(icadv)
        ic = icadv[ind]
        vflux = 0.0
        dvol  = 0.0
        for is in 1:length(g.iscell[ic])
            iface = g.iscell[ic][is]
            vf = vg.fface[iface]
            vp = vg.volpol[iface]
            if g.icface[iface,1] != ic
                vf = -vf; vp = -vp
            end
            vflux += vf; dvol += vp
        end
        # Implicit form with divergence correction
        vg.fractg[ic] = (fract0[ic]*(1.0+0.5*dvol/vg.vcell[ic]) -
                          vflux/vg.vcell[ic]) /
                         (1.0-0.5*dvol/vg.vcell[ic])

        if vg.fractg[ic] < 0
            funder = max(funder, abs(vg.fractg[ic])*vg.vcell[ic])
        end
        if vg.fractg[ic] > 1
            fover = max(fover, (vg.fractg[ic]-1.0)*vg.vcell[ic])
        end

        vg.fractg[ic] ≤ tolfr && (vg.fractg[ic] = 0.0)
        vg.fractg[ic] ≥ (1-tolfr) && (vg.fractg[ic] = 1.0)
    end

    vg.ebound[1] = max(vg.ebound[1], max(funder, fover))
    vg.ebound[2] += max(funder, fover)
    return nothing
end
