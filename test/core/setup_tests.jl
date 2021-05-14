using CryoGrid
using Test
using InteractiveUtils: @which

include("../types.jl")

@testset "Setup" begin
    @testset "3-layer" begin
        grid = Grid(Vector(0.0:10.0:1000.0))
        strat = Stratigraphy(
            -1.0u"m" => Top(TestBoundary()),
            0.0u"m" => Ground(:testground, TestGroundLayer(), TestGroundProcess()),
            1000.0u"m" => Bottom(TestBoundary())
        )
        function checkfields()
            # for each non-error test case, we need to check that all layers are present in uproto and state
            @test hasproperty(setup.uproto,:top)
            @test hasproperty(setup.uproto,:testground)
            @test hasproperty(setup.uproto,:bottom)
            @test hasproperty(setup.cache,:top)
            @test hasproperty(setup.cache,:testground)
            @test hasproperty(setup.cache,:bottom)
            @test hasproperty(setup.meta,:top)
            @test hasproperty(setup.meta,:testground)
            @test hasproperty(setup.meta,:bottom)
        end
        # case: no variables defined
        @test_throws AssertionError setup = CryoGridSetup(strat,grid)
        # case no prognostic variables defined
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
        )
        @test_throws AssertionError setup = CryoGridSetup(strat,grid)
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Prognostic(:x, Float"J", OnGrid(Cells)),
            Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
        )
        setup = CryoGridSetup(strat,grid)
        checkfields()
        @test hasproperty(setup.uproto.testground,:x)
        @test hasproperty(setup.cache.testground,:k)
        # do not allow prognostic layer variables
        CryoGrid.variables(::TestGroundLayer) = (Prognostic(:w,Float"kg",OnGrid(Cells)),)
        @test_throws AssertionError setup = CryoGridSetup(strat,grid)
        # test inclusion of layer variables
        CryoGrid.variables(::TestGroundLayer) = (Diagnostic(:w,Float"kg",OnGrid(Cells)),)
        setup = CryoGridSetup(strat,grid)
        checkfields()
        @test hasproperty(setup.cache.testground,:w)
        # warn when variable declared as prognostic and diagnostic
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Prognostic(:x, Float"J", OnGrid(Cells)),
            Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
            Prognostic(:w, Float"kg", OnGrid(Cells)),
        )
        @test_logs (:warn,r".*declared as both prognostic/algebraic and diagnostic.*") setup = CryoGridSetup(strat,grid)
        checkfields()
        @test hasproperty(setup.uproto.testground,:x)
        @test hasproperty(setup.uproto.testground,:w)
        # allow non-conflicting duplicates
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Prognostic(:x, Float"J", OnGrid(Cells)),
            Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
            Diagnostic(:w, Float"kg", OnGrid(Cells)),
        )
        setup = CryoGridSetup(strat,grid)
        checkfields()
        @test hasproperty(setup.cache.testground,:w)
        # error when conflicting variables declared
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Prognostic(:x, Float"J", OnGrid(Cells)),
            Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
            Diagnostic(:w, Float"kg/m", OnGrid(Edges)),
        )
        @test_throws AssertionError setup = CryoGridSetup(strat,grid)
        # test scalar variables
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Prognostic(:x, Float"J", OnGrid(Cells)),
            Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
            Diagnostic(:a, Float64, Scalar),
        )
        setup = CryoGridSetup(strat,grid)
        checkfields()
        @test hasproperty(setup.cache.testground,:a)
        setup.cache.testground.a.cache.du[1] = 2.0
        @test setup.cache.testground.a.cache.du[1] == 2.0
        # clean-up method definitions (necessary for re-running test set)
        Base.delete_method(@which CryoGrid.variables(TestGroundLayer(),TestGroundProcess()))
        Base.delete_method(@which CryoGrid.variables(TestGroundLayer()))
    end
    @testset "4-layer" begin
        grid = Grid(Vector(0.0:10.0:1000.0))
        strat = Stratigraphy(
            -1.0u"m" => Top(TestBoundary()), (
                0.0u"m" => Ground(:testground1, TestGroundLayer(), TestGroundProcess()),
                100.0u"m" => Ground(:testground2, TestGroundLayer(), TestGroundProcess()),
            ),
            1000.0u"m" => Bottom(TestBoundary())
        )
        CryoGrid.variables(::TestGroundLayer) = (Diagnostic(:w,Float"kg",OnGrid(Cells)),)
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Prognostic(:x, Float"J", OnGrid(Cells)),
            Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
        )
        setup = CryoGridSetup(strat,grid)
        # check for layers
        @test hasproperty(setup.uproto,:top)
        @test hasproperty(setup.uproto,:testground1)
        @test hasproperty(setup.uproto,:testground2)
        @test hasproperty(setup.uproto,:bottom)
        @test hasproperty(setup.cache,:top)
        @test hasproperty(setup.cache,:testground1)
        @test hasproperty(setup.cache,:testground2)
        @test hasproperty(setup.cache,:bottom)
        # check for variables
        @test hasproperty(setup.uproto.testground1, :x)
        @test hasproperty(setup.cache.testground1, :k)
        @test hasproperty(setup.cache.testground1, :w)
        @test hasproperty(setup.uproto.testground2, :x)
        @test hasproperty(setup.cache.testground2, :k)
        @test hasproperty(setup.cache.testground2, :w)
        # clean-up method definitions (necessary for re-running test set)
        Base.delete_method(@which CryoGrid.variables(TestGroundLayer(),TestGroundProcess()))
        Base.delete_method(@which CryoGrid.variables(TestGroundLayer()))
    end
end
