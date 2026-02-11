# ---------------------------------------------------------------------------
#  TESTREC – Reconstruction test driver
#  Julia translation of gvof/testrec.f
# ---------------------------------------------------------------------------
#  Test program to compute the reconstruction errors using different
#  reconstruction methods and grid types.
#
#  Reconstruction test cases:
#    ICASE=101: sphere of radius 0.325 centred at (0.525,0.464,0.516)
#    ICASE=102: sphere of radius 0.4 with a centred spherical hole of
#               radius 0.2, centred at (0.525,0.464,0.516)
#    ICASE=103: torus centred at (0.525,0.464,0.516) with r1=0.1, r2=0.2
#
#  Usage:
#    julia --project=. examples/testrec.jl [gvofvardef_path]
# ---------------------------------------------------------------------------

using gVOF
using ISOAP: isovtkgrid
using Printf

function main(vardeffile::String = "gvofvardef")

    tstart = time()

    # ── Read parameters ──────────────────────────────────────────────────
    params = vardef(vardeffile)
    icase = params.icase

    # ── Construct the grid ───────────────────────────────────────────────
    vg = vofgrid(params)

    # ── Print the grid ───────────────────────────────────────────────────
    printgrid(vg)

    # ── Volume fraction cell initialisation ──────────────────────────────
    func_init = func3dinic(icase)
    initfgrid!(vg, func_init; nc = params.nc, tolfr = params.tolfr)

    # ── Volume fraction node interpolation and tagging ───────────────────
    taggrid!(vg; tolfr = params.tolfr)

    tgrid = time() - tstart

    # ── Interface reconstruction ─────────────────────────────────────────
    tstart = time()

    reconstruct!(vg, params)

    trec = time() - tstart

    # ── Print results ────────────────────────────────────────────────────
    tstart = time()

    params.iprintplic  && printplic(vg, 0)
    if params.iprintiso
        isovtkgrid(vg.grid, vg.frnod, 0.5; ifile = 0)
    end
    params.iprintvoxel && printvoxel(vg, 0)

    tvis = time() - tstart

    # ── Initialisation and reconstruction error ──────────────────────────
    vexact = funcexactvol(icase)
    erec, erinit = recerr(vg, func_init; nc = params.nc, vexact = vexact)

    # ── Report ───────────────────────────────────────────────────────────
    ncellint = length(vg.icint)

    println("-------------------------------------------------")
    println("|            INITIALIZATION ERROR:               |")
    println("-------------------------------------------------")
    @printf("  %.15e\n", erinit)
    println("-------------------------------------------------")
    println("|            RECONSTRUCTION ERROR:               |")
    println("-------------------------------------------------")
    @printf("  %.15e\n", erec)
    println("-------------------------------------------------")
    println("|              EXECUTION TIMES:                  |")
    println("-------------------------------------------------")
    @printf("RECONSTRUCTION = %.4f secs\n", trec)
    @printf("VOF INITIALIZATION AND GRIDDING = %.4f secs\n", tgrid)
    @printf("VISUALIZATION = %.4f secs\n", tvis)
    println("-------------------------------------------------")
    println("NUMBER OF INTERFACIAL CELLS = ", ncellint)
    println("-------------------------------------------------")

    # ── Write log file ───────────────────────────────────────────────────
    open("testrec.logout", "w") do fp
        println(fp, "-------------------------------------------------")
        println(fp, "|            INITIALIZATION ERROR:               |")
        println(fp, "-------------------------------------------------")
        @printf(fp, "  %.15e\n", erinit)
        println(fp, "-------------------------------------------------")
        println(fp, "|            RECONSTRUCTION ERROR:               |")
        println(fp, "-------------------------------------------------")
        @printf(fp, "  %.15e\n", erec)
        println(fp, "-------------------------------------------------")
        println(fp, "|              EXECUTION TIMES:                  |")
        println(fp, "-------------------------------------------------")
        @printf(fp, "RECONSTRUCTION = %.4f secs\n", trec)
        @printf(fp, "VOF INITIALIZATION AND GRIDDING = %.4f secs\n", tgrid)
        println(fp, "-------------------------------------------------")
        println(fp, "NUMBER OF INTERFACIAL CELLS = ", ncellint)
        println(fp, "-------------------------------------------------")
    end

    return (; erinit, erec, ncellint)
end

# ── Run ──────────────────────────────────────────────────────────────────────
if abspath(PROGRAM_FILE) == @__FILE__
    vardeffile = length(ARGS) >= 1 ? ARGS[1] : "gvofvardef"
    main(vardeffile)
end
