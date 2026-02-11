# ---------------------------------------------------------------------------
#  TESTVOF – Advection test driver
#  Julia translation of gvof/testvof.f
# ---------------------------------------------------------------------------
#  Test program to solve the interface dynamics of an initialised material
#  volume body under different prescribed velocity fields and grid types.
#
#  Test cases (advection):
#    ICASE=1: simple translation
#    ICASE=2: XY rotation
#    ICASE=3: 3D deformation  [Enright et al., J. Comput. Phys. 183 (2002)]
#    ICASE=4: vortex in a box XZ [Bell et al., J. Comput. Phys. 85 (1989)]
#    ICASE=5: 3D shearing [Liovic et al., Comput. Fluids 35 (2006)]
#
#  Usage:
#    julia --project=. examples/testvof.jl [gvofvardef_path]
# ---------------------------------------------------------------------------

using gVOF
using ISOAP: isovtkgrid
using Printf

function main(vardeffile::String = "gvofvardef")

    ebound = [0.0, 0.0]   # [max, sum] boundedness errors
    tadv  = 0.0
    trec  = 0.0
    tvis  = 0.0
    tgrid = 0.0

    tstart = time()

    # ── Read parameters ──────────────────────────────────────────────────
    params = vardef(vardeffile)
    icase   = params.icase
    irec    = params.irec
    iadv    = params.iadv
    igrid   = params.igrid
    tcaden  = params.tcaden
    timeend = params.timeend
    tolfr   = params.tolfr
    cfl     = params.cfl
    idt     = params.idt

    toltime = timeend * 1e-16
    icocaden = 1

    # ── Construct the grid ───────────────────────────────────────────────
    vg = vofgrid(params)

    # ── Print the grid ───────────────────────────────────────────────────
    printgrid(vg)

    ncell = vg.grid.ncell

    # ── Allocate advection working arrays ────────────────────────────────
    # fracte stores the exact (end-state) volume fractions
    fracte = zeros(Float64, ncell)

    # ── Volume fraction cell initialisation ──────────────────────────────
    func_init = func3dinic(icase)
    initfgrid!(vg, func_init; nc = params.nc, tolfr = tolfr)

    # ── Total liquid volume ──────────────────────────────────────────────
    vt = sum(vg.fractg[ic] * vg.vcell[ic] for ic in 1:ncell)

    # ── Volume fraction node interpolation and tagging ───────────────────
    taggrid!(vg; tolfr = tolfr)

    tfinish = time()
    tgrid += tfinish - tstart

    # ── Initial interface reconstruction ─────────────────────────────────
    tstart = time()
    reconstruct!(vg, params)
    tfinish = time()
    trec += tfinish - tstart

    # ── Print initial state ──────────────────────────────────────────────
    tstart = time()
    params.iprintplic  && printplic(vg, 0)
    if params.iprintiso
        isovtkgrid(vg.grid, vg.frnod, 0.5; ifile = 0)
    end
    params.iprintvoxel && printvoxel(vg, 0)
    tfinish = time()
    tvis += tfinish - tstart

    # ── Temporal loop ────────────────────────────────────────────────────
    istep  = 0
    timet0 = 0.0
    timet  = 0.0
    dt     = idt == 1 ? 0.0 : params.dt

    # Initial velocity & time step
    velgrid!(vg, icase, timet0, timet)
    if idt == 1
        dt = compdt!(vg, cfl, dt, timeend, timet)
    end

    avgncellint = 0.0

    while true
        istep += 1
        timet0 = timet
        timet  = timet + dt
        tstep  = 0.0

        velgrid!(vg, icase, timet0, timet)

        @printf("%d  Time: %.8e  dt: %.8e\n", istep, timet, dt)

        # ── Face fluxes & advection ──────────────────────────────────
        tstart = time()
        faceflux!(vg; dt = dt, iadv = iadv, igrid = igrid, tolfr = tolfr)
        vofadv!(vg; tolfr = tolfr)

        # Update boundedness errors from vg.ebound
        ebound[1] = max(ebound[1], vg.ebound[1])
        ebound[2] += vg.ebound[2]

        # ── Re-tag ──────────────────────────────────────────────────
        taggrid!(vg; tolfr = tolfr)
        avgncellint += length(vg.icint)

        tfinish = time()
        tstep  += tfinish - tstart
        tadv   += tfinish - tstart

        # ── Re-reconstruct ──────────────────────────────────────────
        tstart = time()
        reconstruct!(vg, params)
        tfinish = time()
        tstep += tfinish - tstart
        trec  += tfinish - tstart

        # ── Periodic output ─────────────────────────────────────────
        if timet > icocaden * tcaden || timet == timeend
            tstart = time()
            params.iprintplic  && printplic(vg, icocaden)
            if params.iprintiso
                isovtkgrid(vg.grid, vg.frnod, 0.5; ifile = icocaden)
            end
            params.iprintvoxel && printvoxel(vg, icocaden)
            icocaden += 1
            tfinish = time()
            tvis += tfinish - tstart
        end

        # ── Check termination ───────────────────────────────────────
        timet > (timeend - toltime) && break

        # ── Next time step ──────────────────────────────────────────
        if idt == 1
            dt = compdt!(vg, cfl, dt, timeend, timet)
        end

        @printf("CPU-TIME PER STEP: %.4f secs\n", tstep)
    end

    # ── Exact volume fraction at the end of the test ─────────────────────
    func_end = func3dendc(icase)
    # Use a temporary VOFGrid-like approach: store exact fractions in fracte
    # by calling initfgrid! with fracte as target
    # We need to use the same grid, so we temporarily swap fractg
    fractg_save = copy(vg.fractg)
    initfgrid!(vg, func_end; nc = params.nc, tolfr = tolfr)
    fracte .= vg.fractg
    vg.fractg .= fractg_save

    # ── Error computation ────────────────────────────────────────────────
    evf = sum(abs(vg.fractg[ic] - fracte[ic]) * vg.vcell[ic] for ic in 1:ncell)
    vtf = sum(vg.fractg[ic] * vg.vcell[ic] for ic in 1:ncell)

    # ── Report ───────────────────────────────────────────────────────────
    println("-------------------------------------------------")
    println("|                ERROR NORMS:                    |")
    println("-------------------------------------------------")
    @printf("GEOMETRIC ERROR = %.15e\n", evf)
    @printf("RELATIVE GEOMETRIC ERROR = %.15e\n", evf / vt)
    @printf("ABSOLUTE CHANGE OF FLUID VOLUME = %.15e\n", abs(vt - vtf))
    @printf("MAXIMUM BOUNDEDNESS VOLUME ERROR = %.15e\n", ebound[1])
    @printf("AVERAGE BOUNDEDNESS VOLUME ERROR = %.15e\n", ebound[2] / istep)
    println("-------------------------------------------------")
    println("|              EXECUTION TIMES:                  |")
    println("-------------------------------------------------")
    @printf("ADVECTION = %.4f secs\n", tadv)
    @printf("RECONSTRUCTION = %.4f secs\n", trec)
    @printf("VOF INITIALIZATION AND GRIDDING = %.4f secs\n", tgrid)
    @printf("VISUALIZATION = %.4f secs\n", tvis)
    println("-------------------------------------------------")
    println("|        AVERAGE VALUES PER TIME STEP:           |")
    println("-------------------------------------------------")
    @printf("ADVECTION AND RECONSTRUCTION EXECUTION TIME = %.4f secs\n",
            (tadv + trec) / istep)
    @printf("NUMBER OF INTERFACIAL CELLS = %d\n",
            round(Int, avgncellint / istep))
    println("-------------------------------------------------")

    # ── Write log file ───────────────────────────────────────────────────
    open("testvof.logout", "w") do fp
        println(fp, "-------------------------------------------------")
        println(fp, "|                ERROR NORMS:                    |")
        println(fp, "-------------------------------------------------")
        @printf(fp, "GEOMETRIC ERROR = %.15e\n", evf)
        @printf(fp, "RELATIVE GEOMETRIC ERROR = %.15e\n", evf / vt)
        @printf(fp, "ABSOLUTE CHANGE OF FLUID VOLUME = %.15e\n", abs(vt - vtf))
        @printf(fp, "MAXIMUM BOUNDEDNESS VOLUME ERROR = %.15e\n", ebound[1])
        @printf(fp, "AVERAGE BOUNDEDNESS VOLUME ERROR = %.15e\n", ebound[2] / istep)
        println(fp, "-------------------------------------------------")
        println(fp, "|              EXECUTION TIMES:                  |")
        println(fp, "-------------------------------------------------")
        @printf(fp, "ADVECTION = %.4f secs\n", tadv)
        @printf(fp, "RECONSTRUCTION = %.4f secs\n", trec)
        @printf(fp, "VOF INITIALIZATION AND GRIDDING = %.4f secs\n", tgrid)
        @printf(fp, "VISUALIZATION = %.4f secs\n", tvis)
        println(fp, "-------------------------------------------------")
        println(fp, "|        AVERAGE VALUES PER TIME STEP:           |")
        println(fp, "-------------------------------------------------")
        @printf(fp, "ADVECTION AND RECONSTRUCTION EXECUTION TIME = %.4f secs\n",
                (tadv + trec) / istep)
        @printf(fp, "NUMBER OF INTERFACIAL CELLS = %d\n",
                round(Int, avgncellint / istep))
        println(fp, "-------------------------------------------------")
    end

    return (; evf, vtf, vt, ebound, istep)
end


# ── Run ──────────────────────────────────────────────────────────────────────
if abspath(PROGRAM_FILE) == @__FILE__
    vardeffile = length(ARGS) >= 1 ? ARGS[1] : "gvofvardef"
    main(vardeffile)
end
