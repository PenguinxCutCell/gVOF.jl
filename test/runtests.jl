using Test
using gVOF

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

    @testset "Reconstruction methods dispatch (irec=1:6)" begin
        for irec in 1:6
            p = VOFParams(nx=6, ny=6, nz=6, irec=irec, niter=2, tolir=1e-2)
            vg = vofgrid(p)
            initfgrid!(vg, func3dinic(101); nc=4, tolfr=p.tolfr)
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
