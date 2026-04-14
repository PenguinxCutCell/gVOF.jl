using Test
using gVOF

function _line_fraction(vg, ic::Int, nx::Float64, ny::Float64, c::Float64)
    poly = gVOF.defcell(vg, ic)
    gVOF.inte3d!(poly, c, nx, ny, 0.0)
    return gVOF.toolv3d(poly) / vg.vcell[ic]
end

function _assign_line_fractions!(vg, nx::Float64, ny::Float64, c::Float64)
    for ic in eachindex(vg.fractg)
        vg.fractg[ic] = _line_fraction(vg, ic, nx, ny, c)
    end
    return nothing
end

@testset "gVOF.jl" begin

    @testset "Types" begin
        p = VOFParams()
        @test p.irec == 1
        @test p.nx == 32
        @test p.xlast == 1.0
    end

    @testset "Grid construction (4×4×4)" begin
        p = VOFParams(nx=4, ny=4, nz=4)
        vg = vofgrid(p)
        @test vg.grid.ncell == 64
        @test vg.grid.npoint == 125
        @test all(vg.vcell .> 0)
        @test isapprox(sum(vg.vcell), 1.0; atol=1e-10)
        @test all(length.(vg.ineigb) .> 0)
    end

    @testset "Test case shapes" begin
        f = func3dinic(1)
        @test f(0.5, 0.5, 0.5) > 0      # inside sphere
        @test f(0.0, 0.0, 0.0) < 0      # outside sphere

        f3 = func3dinic(3)
        @test f3(0.35, 0.35, 0.35) > 0

        vex = funcexactvol(1)
        @test isapprox(vex, 4π*0.25^3/3; rtol=1e-10)
    end

    @testset "Initialisation + Tagging (8×8×8)" begin
        p = VOFParams(nx=8, ny=8, nz=8)
        vg = vofgrid(p)
        f = func3dinic(1)
        vft = initfgrid!(vg, f; nc=p.nc, tolfr=p.tolfr)
        @test vft > 0
        @test any(vg.fractg .> 0) && any(vg.fractg .< 1)

        taggrid!(vg; tolfr=p.tolfr)
        @test length(vg.icint) > 0
        @test length(vg.icadv) > 0
        @test length(vg.isflu) > 0
    end

    @testset "LSGIR reconstruction (8×8×8)" begin
        p = VOFParams(nx=8, ny=8, nz=8, irec=4)
        vg = vofgrid(p)
        initfgrid!(vg, func3dinic(1); nc=p.nc, tolfr=p.tolfr)
        taggrid!(vg; tolfr=p.tolfr)
        lsgir!(vg; pn=p.pn, igrid=p.igrid)
        # Check that normals are set for interfacial cells
        for ic in vg.icint
            d = sqrt(vg.xnormg[ic]^2 + vg.ynormg[ic]^2 + vg.znormg[ic]^2)
            @test isapprox(d, 1.0; atol=1e-10)
        end
    end

    @testset "Reconstruction methods dispatch (irec=1:8)" begin
        for irec in 1:8
            if irec <= 6
                p = VOFParams(nx=6, ny=6, nz=6, irec=irec, niter=2, tolir=1e-2)
                vg = vofgrid(p)
                initfgrid!(vg, func3dinic(101); nc=4, tolfr=p.tolfr)
            else
                p = VOFParams(nx=8, ny=8, nz=1, irec=irec, niter=32, tolir=1e-10)
                vg = vofgrid(p)
                nx = cos(0.37 * π)
                ny = sin(0.37 * π)
                _assign_line_fractions!(vg, nx, ny, -0.08)
            end
            taggrid!(vg; tolfr=p.tolfr)
            @test length(vg.icint) > 0

            reconstruct!(vg, p)

            for ic in vg.icint
                d = sqrt(vg.xnormg[ic]^2 + vg.ynormg[ic]^2 + vg.znormg[ic]^2)
                @test isfinite(vg.rholig[ic])
                @test isfinite(d)
                @test isapprox(d, 1.0; atol=1e-8)
            end
        end
    end

    @testset "ELVIRA and LVIRA exact line reproduction (2D Cartesian)" begin
        cases = (
            (0.17 * π, -0.08),
            (0.41 * π, 0.05),
            (0.63 * π, -0.12),
        )

        for method in (7, 8)
            for (angle, offset) in cases
                nx_exact = cos(angle)
                ny_exact = sin(angle)
                p = VOFParams(nx=16, ny=16, nz=1, irec=method, niter=64, tolir=1e-12)
                vg = vofgrid(p)
                _assign_line_fractions!(vg, nx_exact, ny_exact, offset)
                taggrid!(vg; tolfr=p.tolfr)
                reconstruct!(vg, p)

                cell_of_ij, ij_of_cell = gVOF._build_cartesian_lookup(vg)
                work = gVOF._get_recon_work(vg)

                for ic in vg.icint
                    i, j = ij_of_cell[ic]
                    stencil = gVOF._full_3x3_stencil(cell_of_ij, i, j)
                    stencil === nothing && continue

                    stencil_ids = collect(stencil)
                    objective = gVOF._constrained_stencil_objective!(work, vg, stencil_ids, ic,
                                                                     vg.xnormg[ic], vg.ynormg[ic], vg.znormg[ic],
                                                                     vg.rholig[ic])
                    center_fraction = gVOF._fraction_for_cell_with_plane!(work, vg, ic,
                                                                         vg.xnormg[ic], vg.ynormg[ic], vg.znormg[ic],
                                                                         vg.rholig[ic])
                    dot_err = min(sqrt((vg.xnormg[ic] - nx_exact)^2 + (vg.ynormg[ic] - ny_exact)^2),
                                  sqrt((vg.xnormg[ic] + nx_exact)^2 + (vg.ynormg[ic] + ny_exact)^2))

                    @test isapprox(center_fraction, vg.fractg[ic]; atol=1e-12, rtol=1e-12)
                    @test objective ≤ 1e-10
                    @test dot_err ≤ 1e-7
                    @test isapprox(vg.znormg[ic], 0.0; atol=1e-12)
                end
            end
        end
    end

    @testset "ELVIRA/LVIRA boundary fallback" begin
        for method in (7, 8)
            p = VOFParams(nx=16, ny=16, nz=1, irec=method, niter=64, tolir=1e-12)
            vg = vofgrid(p)
            nx_exact = cos(0.08 * π)
            ny_exact = sin(0.08 * π)
            _assign_line_fractions!(vg, nx_exact, ny_exact, 0.02)
            taggrid!(vg; tolfr=p.tolfr)
            reconstruct!(vg, p)

            @test all(ic -> isfinite(vg.xnormg[ic]) && isfinite(vg.ynormg[ic]) &&
                           isfinite(vg.znormg[ic]) && isfinite(vg.rholig[ic]), vg.icint)
            @test all(ic -> isapprox(sqrt(vg.xnormg[ic]^2 + vg.ynormg[ic]^2 + vg.znormg[ic]^2), 1.0; atol=1e-8),
                      vg.icint)
        end
    end

    @testset "ELVIRA/LVIRA unsupported grid" begin
        for method in (7, 8)
            p = VOFParams(nx=8, ny=8, nz=2, irec=method, niter=16, tolir=1e-10)
            vg = vofgrid(p)
            initfgrid!(vg, func3dinic(1); nc=p.nc, tolfr=p.tolfr)
            taggrid!(vg; tolfr=p.tolfr)
            err = try
                reconstruct!(vg, p)
                nothing
            catch e
                e
            end
            @test err !== nothing
            @test occursin("ELVIRA/LVIRA currently support only 2D uniform Cartesian grids with nz == 1",
                           sprint(showerror, err))
        end
    end

    @testset "ELVIRA circle convergence (2D Cartesian)" begin
        errs = Float64[]
        circle = (x, y, z) -> 0.22^2 - ((x - 0.5)^2 + (y - 0.5)^2)
        for n in (32, 64, 128)
            p = VOFParams(nx=n, ny=n, nz=1, irec=7, nc=6, tolfr=1e-12)
            vg = vofgrid(p)
            initfgrid!(vg, circle; nc=p.nc, tolfr=p.tolfr)
            taggrid!(vg; tolfr=p.tolfr)
            reconstruct!(vg, p)
            err = recerr(vg, circle; nc=p.nc, vexact=π * 0.22^2)
            push!(errs, err.erec)
        end
        @test all(isfinite, errs)
        @test errs[2] < errs[1]
        @test errs[3] < errs[2]
    end

    @testset "Advection one-step sanity (translation case)" begin
        p = VOFParams(; icase=1, irec=4, iadv=1, nx=6, ny=6, nz=6,
                      idt=0, dt=0.01, tolfr=1e-12)
        vg = vofgrid(p)
        initfgrid!(vg, func3dinic(p.icase); nc=4, tolfr=p.tolfr)
        taggrid!(vg; tolfr=p.tolfr)
        reconstruct!(vg, p)

        v0 = sum(vg.fractg .* vg.vcell)

        velgrid!(vg, p.icase, 0.0, p.dt)
        faceflux!(vg; dt=p.dt, iadv=p.iadv, igrid=p.igrid, tolfr=p.tolfr)
        vofadv!(vg; tolfr=p.tolfr)
        taggrid!(vg; tolfr=p.tolfr)
        reconstruct!(vg, p)

        v1 = sum(vg.fractg .* vg.vcell)
        @test all(0.0 .<= vg.fractg .<= 1.0)
        @test abs(v1 - v0) ≤ 2e-2
        @test vg.ebound[1] ≤ 1e-6
    end

    @testset "gvofvardef-style geometric error (no file I/O)" begin
        # Matches examples/gvofvardef settings, with reduced grid size for test speed.
        p = VOFParams(
            icase=3, irec=1, pn=1.5, iw=1, niter=4, tolir=1e-3,
            iadv=1, tcaden=0.1, iprintplic=true, iprintiso=false, iprintvoxel=false,
            nc=4, tolfr=1e-12, idt=1, dt=1e-5, cfl=0.5, igrid=1,
            nx=8, ny=8, nz=8, xlast=1.0, ylast=1.0, zlast=1.0, timeend=3.0,
        )
        vg = vofgrid(p)
        f = func3dinic(p.icase)

        initfgrid!(vg, f, p.nc, p.tolfr)
        taggrid!(vg, p.tolfr)
        reconstruct!(vg, p)

        err = recerr(vg, f; nc=p.nc, vexact=funcexactvol(p.icase))
        @test length(vg.icint) > 0
        @test isfinite(err.erec) && err.erec ≥ 0.0
        @test isfinite(err.erinit) && err.erinit ≥ 0.0
        @test err.erec ≤ 5e-3
        @test err.erinit ≤ 5e-2
    end

end
