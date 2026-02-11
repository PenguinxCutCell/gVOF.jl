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

end
