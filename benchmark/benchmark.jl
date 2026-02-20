#!/usr/bin/env julia

using gVOF
using Printf
using Statistics
using Profile

"""
Simple benchmarking script for gVOF kernels.

Usage examples:
  julia --project=. benchmark/benchmark.jl
  julia --project=. benchmark/benchmark.jl nx=20 reps=4 icase=3
  julia --project=. benchmark/benchmark.jl nx=24 reps=3 do_profile=true
"""

function parse_args(args)
    cfg = Dict(
        "nx" => 12,
        "ny" => 12,
        "nz" => 12,
        "icase" => 3,
        "irec" => 4,
        "iadv" => 1,
        "reps" => 3,
        "nc" => 6,
        "dt" => 0.01,
        "tolfr" => 1e-12,
        "do_profile" => false,
    )

    for a in args
        occursin('=', a) || continue
        k, v = split(a, '='; limit=2)
        if !haskey(cfg, k)
            @warn "Ignoring unknown option" k
            continue
        end
        if cfg[k] isa Bool
            cfg[k] = lowercase(v) in ("1", "true", "yes", "on")
        elseif cfg[k] isa Int
            cfg[k] = parse(Int, v)
        elseif cfg[k] isa Float64
            cfg[k] = parse(Float64, v)
        end
    end
    return cfg
end

function run_repeated(f; reps::Int)
    times = Float64[]
    bytes = Int[]
    gctimes = Float64[]

    for _ in 1:reps
        GC.gc()
        s = @timed f()
        push!(times, s.time)
        push!(bytes, s.bytes)
        push!(gctimes, s.gctime)
    end
    return (
        t_mean = mean(times),
        t_min = minimum(times),
        t_max = maximum(times),
        b_mean = mean(bytes),
        gc_mean = mean(gctimes),
    )
end

function print_row(name::AbstractString, m)
    @printf("%-34s %10.4f %10.4f %10.4f %12.1f %10.4f\n",
            name, 1e3*m.t_mean, 1e3*m.t_min, 1e3*m.t_max, m.b_mean/1024^2, 1e3*m.gc_mean)
end

function print_header(title::AbstractString)
    println()
    println(title)
    println("-"^95)
    @printf("%-34s %10s %10s %10s %12s %10s\n",
            "Kernel", "mean ms", "min ms", "max ms", "mean MiB", "GC ms")
    println("-"^95)
end

function fresh_state(p::VOFParams; nc::Int, tolfr::Float64)
    vg = vofgrid(p)
    initfgrid!(vg, func3dinic(p.icase), nc, tolfr)
    taggrid!(vg, tolfr)
    reconstruct!(vg, p)
    return vg
end

function benchmark_pipeline(cfg)
    p = VOFParams(
        icase = cfg["icase"],
        irec = cfg["irec"],
        iadv = cfg["iadv"],
        nx = cfg["nx"], ny = cfg["ny"], nz = cfg["nz"],
        idt = 0, dt = cfg["dt"], tolfr = cfg["tolfr"],
    )
    reps = cfg["reps"]
    nc = cfg["nc"]
    tolfr = cfg["tolfr"]

    print_header("Pipeline breakdown (single configuration)")
    @printf("Config: icase=%d irec=%d iadv=%d nx=%d ny=%d nz=%d reps=%d nc=%d dt=%.4g\n",
            p.icase, p.irec, p.iadv, p.nx, p.ny, p.nz, reps, nc, p.dt)

    m_grid = run_repeated(() -> vofgrid(p); reps=reps)
    print_row("vofgrid", m_grid)

    m_init_tag = run_repeated(() -> begin
        vg = vofgrid(p)
        initfgrid!(vg, func3dinic(p.icase), nc, tolfr)
        taggrid!(vg, tolfr)
        nothing
    end; reps=reps)
    print_row("initfgrid! + taggrid!", m_init_tag)

    m_reconstruct = run_repeated(() -> begin
        vg = vofgrid(p)
        initfgrid!(vg, func3dinic(p.icase), nc, tolfr)
        taggrid!(vg, tolfr)
        reconstruct!(vg, p)
        nothing
    end; reps=reps)
    print_row("reconstruct!", m_reconstruct)

    m_adv_step = run_repeated(() -> begin
        vg = fresh_state(p; nc=nc, tolfr=tolfr)
        velgrid!(vg, p.icase, 0.0, p.dt)
        faceflux!(vg; dt=p.dt, iadv=p.iadv, igrid=p.igrid, tolfr=tolfr)
        vofadv!(vg, tolfr)
        taggrid!(vg, tolfr)
        reconstruct!(vg, p)
        nothing
    end; reps=reps)
    print_row("one advection step", m_adv_step)
    println("-"^95)
end

function benchmark_reconstruction_sweep(cfg)
    reps = cfg["reps"]
    nc = cfg["nc"]
    tolfr = cfg["tolfr"]

    print_header("Reconstruction sweep (irec = 1:6)")
    for irec in 1:6
        p = VOFParams(
            icase = 101,
            irec = irec,
            nx = cfg["nx"], ny = cfg["ny"], nz = cfg["nz"],
            tolfr = tolfr,
            niter = 2,
            tolir = 1e-2,
        )

        m = run_repeated(() -> begin
            vg = vofgrid(p)
            initfgrid!(vg, func3dinic(101), nc, tolfr)
            taggrid!(vg, tolfr)
            reconstruct!(vg, p)
            nothing
        end; reps=reps)

        print_row("irec=$irec reconstruction", m)
    end
    println("-"^95)
end

function benchmark_advection_sweep(cfg)
    reps = cfg["reps"]
    nc = cfg["nc"]
    tolfr = cfg["tolfr"]

    print_header("Advection sweep (iadv = 1:3, one step)")
    for iadv in 1:3
        p = VOFParams(
            icase = 1,
            irec = 4,
            iadv = iadv,
            nx = cfg["nx"], ny = cfg["ny"], nz = cfg["nz"],
            idt = 0, dt = cfg["dt"], tolfr = tolfr,
        )

        m = run_repeated(() -> begin
            vg = fresh_state(p; nc=nc, tolfr=tolfr)
            velgrid!(vg, p.icase, 0.0, p.dt)
            faceflux!(vg; dt=p.dt, iadv=iadv, igrid=p.igrid, tolfr=tolfr)
            vofadv!(vg, tolfr)
            taggrid!(vg, tolfr)
            reconstruct!(vg, p)
            nothing
        end; reps=reps)

        print_row("iadv=$iadv one-step", m)
    end
    println("-"^95)
end

function run_profile_snapshot(cfg)
    do_profile = cfg["do_profile"]
    do_profile || return

    p = VOFParams(
        icase = cfg["icase"],
        irec = cfg["irec"],
        iadv = cfg["iadv"],
        nx = cfg["nx"], ny = cfg["ny"], nz = cfg["nz"],
        idt = 0, dt = cfg["dt"], tolfr = cfg["tolfr"],
    )
    nc = cfg["nc"]
    tolfr = cfg["tolfr"]

    println()
    println("Profile snapshot (flat, by sample count)")
    println("-"^95)

    vg = fresh_state(p; nc=nc, tolfr=tolfr)
    velgrid!(vg, p.icase, 0.0, p.dt)
    Profile.clear()
    @profile begin
        faceflux!(vg; dt=p.dt, iadv=p.iadv, igrid=p.igrid, tolfr=tolfr)
        vofadv!(vg, tolfr)
        taggrid!(vg, tolfr)
        reconstruct!(vg, p)
    end
    Profile.print(format=:flat, sortedby=:count, maxdepth=12)
    println("-"^95)
end

function main(args=ARGS)
    cfg = parse_args(args)
    benchmark_pipeline(cfg)
    benchmark_reconstruction_sweep(cfg)
    benchmark_advection_sweep(cfg)
    run_profile_snapshot(cfg)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
