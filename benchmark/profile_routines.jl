#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using gVOF
using Profile
using Base.Threads

function _base_params(; nx::Int=12, ny::Int=12, nz::Int=12,
                      icase::Int=3, irec::Int=4, iadv::Int=1, dt::Float64=0.01)
    return VOFParams(
        icase = icase,
        irec = irec,
        iadv = iadv,
        nx = nx,
        ny = ny,
        nz = nz,
        idt = 0,
        dt = dt,
        tolfr = 1e-12,
        niter = 2,
        tolir = 1e-2,
    )
end

function _prepare_vg_for_tagging(p::VOFParams; nc::Int=6)
    vg = vofgrid(p)
    initfgrid!(vg, func3dinic(p.icase), nc, p.tolfr)
    return vg
end

function _prepare_vg_for_reconstruct(p::VOFParams; nc::Int=6)
    vg = _prepare_vg_for_tagging(p; nc=nc)
    taggrid!(vg, p.tolfr)
    return vg
end

function _prepare_vg_for_advection(p::VOFParams; nc::Int=6, timet0::Float64=0.0)
    vg = _prepare_vg_for_reconstruct(p; nc=nc)
    reconstruct!(vg, p)
    velgrid!(vg, p.icase, timet0, timet0 + p.dt)
    return vg
end

function build_case(case::Symbol; nx::Int=12, ny::Int=12, nz::Int=12,
                    icase::Int=3, irec::Int=4, iadv::Int=1, dt::Float64=0.01,
                    nc::Int=6)
    if case == :vofgrid
        p = _base_params(; nx=nx, ny=ny, nz=nz, icase=icase, irec=irec, iadv=iadv, dt=dt)
        vg = vofgrid(p)
        return let vg = vg, p = p
            () -> vofgrid!(vg, p)
        end
    elseif case == :vofgrid_alloc
        p = _base_params(; nx=nx, ny=ny, nz=nz, icase=icase, irec=irec, iadv=iadv, dt=dt)
        return let p = p
            () -> vofgrid(p)
        end
    elseif case == :initfgrid
        p = _base_params(; nx=nx, ny=ny, nz=nz, icase=icase, irec=irec, iadv=iadv, dt=dt)
        vg = vofgrid(p)
        f = func3dinic(p.icase)
        tolfr = p.tolfr
        return let vg = vg, f = f, nc = nc, tolfr = tolfr
            () -> initfgrid!(vg, f, nc, tolfr)
        end
    elseif case == :taggrid
        p = _base_params(; nx=nx, ny=ny, nz=nz, icase=icase, irec=irec, iadv=iadv, dt=dt)
        vg = _prepare_vg_for_tagging(p; nc=nc)
        tolfr = p.tolfr
        return let vg = vg, tolfr = tolfr
            () -> taggrid!(vg, tolfr)
        end
    elseif case == :reconstruct
        p = _base_params(; nx=nx, ny=ny, nz=nz, icase=icase, irec=irec, iadv=iadv, dt=dt)
        vg = _prepare_vg_for_reconstruct(p; nc=nc)
        return let vg = vg, p = p
            () -> reconstruct!(vg, p)
        end
    elseif case == :reconstruct_all
        p0 = _base_params(; nx=nx, ny=ny, nz=nz, icase=101, irec=1, iadv=iadv, dt=dt)
        vg = _prepare_vg_for_reconstruct(p0; nc=nc)
        return let vg = vg, p0 = p0
            () -> begin
                for m in 1:6
                    p2 = VOFParams(
                        icase = p0.icase,
                        irec = m,
                        pn = p0.pn,
                        iw = p0.iw,
                        niter = p0.niter,
                        tolir = p0.tolir,
                        iadv = p0.iadv,
                        tcaden = p0.tcaden,
                        iprintplic = p0.iprintplic,
                        iprintiso = p0.iprintiso,
                        iprintvoxel = p0.iprintvoxel,
                        nc = p0.nc,
                        tolfr = p0.tolfr,
                        idt = p0.idt,
                        dt = p0.dt,
                        cfl = p0.cfl,
                        igrid = p0.igrid,
                        nx = p0.nx,
                        ny = p0.ny,
                        nz = p0.nz,
                        xlast = p0.xlast,
                        ylast = p0.ylast,
                        zlast = p0.zlast,
                        timeend = p0.timeend,
                    )
                    reconstruct!(vg, p2)
                end
                nothing
            end
        end
    elseif case == :faceflux
        p = _base_params(; nx=nx, ny=ny, nz=nz, icase=icase, irec=irec, iadv=iadv, dt=dt)
        vg = _prepare_vg_for_advection(p; nc=nc)
        dt2 = p.dt
        iadv2 = p.iadv
        igrid2 = p.igrid
        tolfr = p.tolfr
        return let vg = vg, dt2 = dt2, iadv2 = iadv2, igrid2 = igrid2, tolfr = tolfr
            () -> faceflux!(vg; dt=dt2, iadv=iadv2, igrid=igrid2, tolfr=tolfr)
        end
    elseif case == :vofadv
        p = _base_params(; nx=nx, ny=ny, nz=nz, icase=icase, irec=irec, iadv=iadv, dt=dt)
        vg = _prepare_vg_for_advection(p; nc=nc)
        faceflux!(vg; dt=p.dt, iadv=p.iadv, igrid=p.igrid, tolfr=p.tolfr)
        fract0 = copy(vg.fractg)
        tolfr = p.tolfr
        return let vg = vg, fract0 = fract0, tolfr = tolfr
            () -> begin
                vg.fractg .= fract0
                vg.ebound .= 0.0
                vofadv!(vg, tolfr)
            end
        end
    elseif case == :advection_step
        p = _base_params(; nx=nx, ny=ny, nz=nz, icase=icase, irec=irec, iadv=iadv, dt=dt)
        vg = _prepare_vg_for_advection(p; nc=nc)
        fract0 = copy(vg.fractg)
        dt2 = p.dt
        iadv2 = p.iadv
        igrid2 = p.igrid
        tolfr = p.tolfr
        icase2 = p.icase
        return let vg = vg, fract0 = fract0, dt2 = dt2, iadv2 = iadv2,
                   igrid2 = igrid2, tolfr = tolfr, icase2 = icase2, p = p
            () -> begin
                vg.fractg .= fract0
                vg.ebound .= 0.0
                velgrid!(vg, icase2, 0.0, dt2)
                faceflux!(vg; dt=dt2, iadv=iadv2, igrid=igrid2, tolfr=tolfr)
                vofadv!(vg, tolfr)
                taggrid!(vg, tolfr)
                reconstruct!(vg, p)
                nothing
            end
        end
    elseif case == :pipeline
        p = _base_params(; nx=nx, ny=ny, nz=nz, icase=icase, irec=irec, iadv=iadv, dt=dt)
        f = func3dinic(p.icase)
        vg = vofgrid(p)
        dt2 = p.dt
        iadv2 = p.iadv
        igrid2 = p.igrid
        tolfr = p.tolfr
        icase2 = p.icase
        return let p = p, vg = vg, f = f, nc = nc, dt2 = dt2, iadv2 = iadv2,
                   igrid2 = igrid2, tolfr = tolfr, icase2 = icase2
            () -> begin
                vofgrid!(vg, p)
                initfgrid!(vg, f, nc, tolfr)
                taggrid!(vg, tolfr)
                reconstruct!(vg, p)
                velgrid!(vg, icase2, 0.0, dt2)
                faceflux!(vg; dt=dt2, iadv=iadv2, igrid=igrid2, tolfr=tolfr)
                vofadv!(vg, tolfr)
                taggrid!(vg, tolfr)
                reconstruct!(vg, p)
                nothing
            end
        end
    else
        error("Unknown case: $(case). Use one of :vofgrid, :vofgrid_alloc, :initfgrid, :taggrid, :reconstruct, :reconstruct_all, :faceflux, :vofadv, :advection_step, :pipeline.")
    end
end

function benchmark_profile(f::F; case::Symbol=:reconstruct, n_eval::Int=100,
                           profile_n::Int=25, nx::Int=12, ny::Int=12, nz::Int=12,
                           icase::Int=3, irec::Int=4, iadv::Int=1, dt::Float64=0.01,
                           nc::Int=6) where {F}
    # Warm-up compilation before measuring runtime/allocations.
    f()

    elapsed = @elapsed begin
        for _ in 1:n_eval
            f()
        end
    end

    total_alloc_bytes = @allocated begin
        for _ in 1:n_eval
            f()
        end
    end

    println("Case: $(case)")
    println("Threads: $(nthreads())")
    println("Runs: $(n_eval)")
    println("Grid: nx=$(nx) ny=$(ny) nz=$(nz)")
    println("Params: icase=$(icase) irec=$(irec) iadv=$(iadv) nc=$(nc) dt=$(dt)")
    println("Total time (s): ", elapsed)
    println("Avg time (us/run): ", elapsed * 1e6 / n_eval)
    println("Total allocated (bytes): ", total_alloc_bytes)
    println("Avg allocated (bytes/run): ", total_alloc_bytes / n_eval)

    Profile.clear()
    @profile begin
        for _ in 1:profile_n
            f()
        end
    end

    println("\nCPU profile (flat, by sample count):")
    Profile.print(format=:flat, sortedby=:count, mincount=5)

    return nothing
end

function benchmark_profile(; case::Symbol=:reconstruct, n_eval::Int=100,
                           profile_n::Int=25, nx::Int=12, ny::Int=12, nz::Int=12,
                           icase::Int=3, irec::Int=4, iadv::Int=1, dt::Float64=0.01,
                           nc::Int=6)
    f = build_case(case; nx=nx, ny=ny, nz=nz, icase=icase, irec=irec, iadv=iadv, dt=dt, nc=nc)
    return benchmark_profile(
        f;
        case=case,
        n_eval=n_eval,
        profile_n=profile_n,
        nx=nx,
        ny=ny,
        nz=nz,
        icase=icase,
        irec=irec,
        iadv=iadv,
        dt=dt,
        nc=nc,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    case      = length(ARGS) >= 1 ? Symbol(ARGS[1]) : :reconstruct
    n_eval    = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 100
    profile_n = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 25
    nx        = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 12
    ny        = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : nx
    nz        = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : nx
    icase     = length(ARGS) >= 7 ? parse(Int, ARGS[7]) : 3
    irec      = length(ARGS) >= 8 ? parse(Int, ARGS[8]) : 4
    iadv      = length(ARGS) >= 9 ? parse(Int, ARGS[9]) : 1
    nc        = length(ARGS) >= 10 ? parse(Int, ARGS[10]) : 6
    dt        = length(ARGS) >= 11 ? parse(Float64, ARGS[11]) : 0.01

    benchmark_profile(; case=case, n_eval=n_eval, profile_n=profile_n,
                      nx=nx, ny=ny, nz=nz, icase=icase, irec=irec, iadv=iadv, nc=nc, dt=dt)
end
