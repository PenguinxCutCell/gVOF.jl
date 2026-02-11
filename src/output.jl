# ---------------------------------------------------------------------------
#  VTK output routines  (PRINTGRID, PRINTPLIC, PRINTVOXEL, CARACTER
#  from gvof.f)
# ---------------------------------------------------------------------------

"""
    _caracter(i::Int) -> String

Convert an integer 0–9999 to a zero-padded 4-character string.
"""
_caracter(i::Int) = @sprintf "%04d" i

# ═══════════════════════════════════════════════════════════════════════════
#  PRINTGRID – write the grid faces as VTK POLYDATA
# ═══════════════════════════════════════════════════════════════════════════
"""
    printgrid(vg::VOFGrid; filename="grid.vtk")

Write all grid faces to a VTK legacy-ASCII POLYDATA file.
"""
function printgrid(vg::VOFGrid; filename::String="grid.vtk")
    g = vg.grid
    npoint = g.npoint
    nface  = g.nface

    open(filename, "w") do io
        println(io, "# vtk DataFile Version 2.0")
        println(io, "Mesh")
        println(io, "ASCII")
        println(io, "DATASET POLYDATA")
        @printf(io, "POINTS %d float\n", npoint)
        for ip in 1:npoint
            @printf(io, "%e %e %e\n", g.vnode[ip,1], g.vnode[ip,2], g.vnode[ip,3])
        end

        # Count total connectivity entries
        ntotal = 0
        for iface in 1:nface
            ntotal += length(g.ipface[iface]) + 1
        end
        @printf(io, "POLYGONS %d %d\n", nface, ntotal)
        for iface in 1:nface
            ipf = g.ipface[iface]
            nip = length(ipf)
            print(io, nip)
            for ip in ipf
                print(io, " ", ip-1)  # 0-based
            end
            println(io)
        end
    end
    @info "Wrote grid to $filename"
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════
#  PRINTPLIC – write reconstructed PLIC interface polygons as VTK POLYDATA
# ═══════════════════════════════════════════════════════════════════════════
"""
    printplic(vg::VOFGrid, ifile::Int; dir=".")

Write the PLIC interface polygons for the current reconstruction to
`plic-<ifile>.vtk` in VTK legacy-ASCII POLYDATA format.
"""
function printplic(vg::VOFGrid, ifile::Int; dir::String=".")
    g     = vg.grid
    icint = vg.icint

    # Collect PLIC polygon vertices and faces
    all_verts = Float64[]
    all_faces = Vector{Int}[]
    voff = 0  # vertex offset

    for ind in eachindex(icint)
        ic = icint[ind]
        poly = defcell(vg, ic)
        xnc = vg.xnormg[ic]; ync = vg.ynormg[ic]; znc = vg.znormg[ic]
        c = vg.rholig[ic]
        nts_before = poly.nts

        icontn, icontp = inte3d!(poly, c, xnc, ync, znc)
        icontp == 0 && continue

        # The newly added face is the PLIC polygon
        if poly.nts > nts_before
            is = poly.nts
            nv = poly.nipv[is]
            local_verts = Int[]
            for iv in 1:nv
                j = poly.ipv[is, iv]
                push!(all_verts, poly.vertp[j,1], poly.vertp[j,2], poly.vertp[j,3])
                push!(local_verts, voff)
                voff += 1
            end
            push!(all_faces, local_verts)
        end
    end

    npt = voff
    nf  = length(all_faces)
    fname = joinpath(dir, "plic-$(_caracter(ifile)).vtk")

    open(fname, "w") do io
        println(io, "# vtk DataFile Version 2.0")
        println(io, "PLIC")
        println(io, "ASCII")
        println(io, "DATASET POLYDATA")
        @printf(io, "POINTS %d float\n", npt)
        for i in 1:3:length(all_verts)
            @printf(io, "%e %e %e\n", all_verts[i], all_verts[i+1], all_verts[i+2])
        end
        ntotal = sum(length(f)+1 for f in all_faces; init=0)
        @printf(io, "POLYGONS %d %d\n", nf, ntotal)
        for f in all_faces
            print(io, length(f))
            for v in f; print(io, " ", v); end
            println(io)
        end
    end
    @info "Wrote PLIC to $fname"
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════
#  PRINTVOXEL – voxelised volume fraction as VTK UNSTRUCTURED_GRID
# ═══════════════════════════════════════════════════════════════════════════
"""
    printvoxel(vg::VOFGrid, ifile::Int; dir=".", nvx=0, nvy=0, nvz=0)

Write a voxel-based volume-fraction field to `voxel-<ifile>.vtk` in
VTK legacy-ASCII UNSTRUCTURED_GRID format.

For a uniform grid (`igrid=1`) the voxel grid matches the computational
grid.  For other grid types a uniform overlay is used.
"""
function printvoxel(vg::VOFGrid, ifile::Int; dir::String=".",
                    nvx::Int=0, nvy::Int=0, nvz::Int=0)
    g = vg.grid

    # Determine voxel grid dimensions from the bounding box
    xmin = minimum(@view vg.boxcell[:,1]); xmax = maximum(@view vg.boxcell[:,2])
    ymin = minimum(@view vg.boxcell[:,3]); ymax = maximum(@view vg.boxcell[:,4])
    zmin = minimum(@view vg.boxcell[:,5]); zmax = maximum(@view vg.boxcell[:,6])

    if nvx == 0
        # Default: use approximate cell count in each direction
        sc = vg.sizecmin
        nvx = max(2, round(Int, (xmax-xmin)/sc[1]))
        nvy = max(2, round(Int, (ymax-ymin)/sc[2]))
        nvz = max(2, round(Int, (zmax-zmin)/sc[3]))
    end

    dx = (xmax-xmin)/nvx; dy = (ymax-ymin)/nvy; dz = (zmax-zmin)/nvz

    # Voxel cell volume fractions – simple nearest-cell approach
    fvox = zeros(nvx, nvy, nvz)
    for ic in 1:g.ncell
        # Find voxels overlapping this cell (using bounding box)
        ix1 = max(1, floor(Int, (vg.boxcell[ic,1]-xmin)/dx) + 1)
        ix2 = min(nvx, floor(Int, (vg.boxcell[ic,2]-xmin)/dx) + 1)
        iy1 = max(1, floor(Int, (vg.boxcell[ic,3]-ymin)/dy) + 1)
        iy2 = min(nvy, floor(Int, (vg.boxcell[ic,4]-ymin)/dy) + 1)
        iz1 = max(1, floor(Int, (vg.boxcell[ic,5]-zmin)/dz) + 1)
        iz2 = min(nvz, floor(Int, (vg.boxcell[ic,6]-zmin)/dz) + 1)
        for iz in iz1:iz2, iy in iy1:iy2, ix in ix1:ix2
            fvox[ix,iy,iz] = vg.fractg[ic]
        end
    end

    # Interpolate to node-centred values
    fnode = zeros(nvx+1, nvy+1, nvz+1)
    for iz in 1:nvz+1, iy in 1:nvy+1, ix in 1:nvx+1
        n = 0; s = 0.0
        for diz in -1:0, diy in -1:0, dix in -1:0
            jx = ix+dix; jy = iy+diy; jz = iz+diz
            (1 ≤ jx ≤ nvx && 1 ≤ jy ≤ nvy && 1 ≤ jz ≤ nvz) || continue
            s += fvox[jx,jy,jz]; n += 1
        end
        fnode[ix,iy,iz] = n > 0 ? s/n : 0.0
    end

    # Write VTK
    npts  = (nvx+1)*(nvy+1)*(nvz+1)
    ncells = nvx*nvy*nvz
    fname = joinpath(dir, "voxel-$(_caracter(ifile)).vtk")

    open(fname, "w") do io
        println(io, "# vtk DataFile Version 2.0")
        println(io, "SCALAR")
        println(io, "ASCII")
        println(io, "DATASET UNSTRUCTURED_GRID")
        @printf(io, "POINTS %d float\n", npts)
        for iz in 0:nvz, iy in 0:nvy, ix in 0:nvx
            @printf(io, "%e %e %e\n", xmin+ix*dx, ymin+iy*dy, zmin+iz*dz)
        end

        @printf(io, "CELLS %d %d\n", ncells, ncells*9)
        for iz in 0:nvz-1, iy in 0:nvy-1, ix in 0:nvx-1
            n0 = iz*(nvy+1)*(nvx+1) + iy*(nvx+1) + ix
            @printf(io, "8 %d %d %d %d %d %d %d %d\n",
                    n0, n0+1,
                    n0+(nvx+1)+1, n0+(nvx+1),
                    n0+(nvy+1)*(nvx+1), n0+(nvy+1)*(nvx+1)+1,
                    n0+(nvy+1)*(nvx+1)+(nvx+1)+1,
                    n0+(nvy+1)*(nvx+1)+(nvx+1))
        end

        @printf(io, "CELL_TYPES %d\n", ncells)
        for _ in 1:ncells
            println(io, "11")
        end

        @printf(io, "POINT_DATA %d\n", npts)
        println(io, "SCALARS vof float 1")
        println(io, "LOOKUP_TABLE default")
        for iz in 1:nvz+1, iy in 1:nvy+1, ix in 1:nvx+1
            @printf(io, "%e\n", fnode[ix,iy,iz])
        end
    end
    @info "Wrote voxel to $fname"
    return nothing
end
