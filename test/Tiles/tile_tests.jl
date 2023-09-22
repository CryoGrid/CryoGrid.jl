using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils
using Test
using InteractiveUtils: @which

include("../types.jl")

@testset "Tile" begin
    @testset "3-layer" begin
        grid = Grid(Vector(0.0:10.0:1000.0)u"m")
        strat = Stratigraphy(
            -1.0u"m" => Top(TestBoundary()),
            0.0u"m" => Named(:testground, TestGroundLayer(TestGroundProcess())),
            1000.0u"m" => Bottom(TestBoundary())
        )
        function checkfields(model)
            # for each non-error test case, we need to check that all layers are present in vars
            @test hasproperty(model.state.vars, :top)
            @test hasproperty(model.state.vars, :testground)
            @test hasproperty(model.state.vars, :bottom)
        end
        try
            # case: no variables defined
            @test_throws AssertionError model = Tile(strat, grid)
            # case: no prognostic variables defined
            CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
                Diagnostic(:k, OnGrid(Edges), u"J/s/m^3"),
            )
            @test_throws AssertionError model = Tile(strat, grid)
            # case: OK
            CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
                Prognostic(:x, OnGrid(Cells), u"J"),
                Diagnostic(:k, OnGrid(Edges), u"J/s/m^3"),
            )
            model = Tile(strat, grid, DummyInitializer{:x}())
            checkfields(model)
            @test hasproperty(model.state.uproto,:x)
            @test hasproperty(model.state.griddiag,:k)
            # no initializers
            model = @test_logs (:warn, r"No initializers provided.*") Tile(strat, grid)
            # variable declared as prognostic and diagnostic
            CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
                Prognostic(:x, OnGrid(Cells), u"J"),
                Diagnostic(:k, OnGrid(Edges), u"J/s/m^3"),
                Diagnostic(:w, OnGrid(Cells), u"kg"),
                Prognostic(:w, OnGrid(Cells), u"kg"),
            )
            model = Tile(strat, grid, DummyInitializer{:x}())
            checkfields(model)
            @test hasproperty(model.state.uproto,:x)
            @test hasproperty(model.state.uproto,:w)
            # allow non-conflicting duplicates
            CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
                Prognostic(:x, OnGrid(Cells), u"J"),
                Diagnostic(:k, OnGrid(Edges), u"J/s/m^3"),
                Diagnostic(:w, OnGrid(Cells), u"kg"),
                Diagnostic(:w, OnGrid(Cells), u"kg"),
            )
            model = Tile(strat, grid, DummyInitializer{:x}())
            checkfields(model)
            @test hasproperty(model.state.griddiag,:w)
            # error when conflicting variables declared
            CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
                Prognostic(:x, OnGrid(Cells), u"J"),
                Diagnostic(:k, OnGrid(Edges), u"J/s/m^3"),
                Diagnostic(:w, OnGrid(Cells), u"kg/m"),
                Diagnostic(:w, OnGrid(Edges), u"kg/m"),
            )
            @test_throws AssertionError model = Tile(strat, grid, DummyInitializer{:x}())
            # test scalar variables and parameters
            CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
                Prognostic(:x, OnGrid(Cells), u"J"),
                Diagnostic(:k, OnGrid(Edges), u"J/s/m^3"),
                Diagnostic(:a, Scalar, NoUnits, Float64),
            )
            model = Tile(strat, grid, DummyInitializer{:x}())
            checkfields(model)
            @test hasproperty(model.state.diag.testground, :a)
            model.state.diag.testground.a.cache.du[1] = 2.0
            @test model.state.diag.testground.a.cache.du[1] == 2.0
            state = getstate(:testground, model, model.state.uproto, model.state.uproto, 0.0)
        finally
            # clean-up method definitions (necessary for re-running test set)
            Base.delete_method(@which CryoGrid.variables(TestGroundLayer(TestGroundProcess()), TestGroundProcess()))
        end
    end
    @testset "4-layer" begin
        grid = Grid(Vector(0.0:10.0:1000.0)u"m")
        strat = Stratigraphy(
            -1.0u"m" => Top(TestBoundary()), (
                0.0u"m" => Named(:testground1, TestGroundLayer(TestGroundProcess())),
                100.0u"m" => Named(:testground2, TestGroundLayer(TestGroundProcess())),
            ),
            1000.0u"m" => Bottom(TestBoundary())
        )
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Diagnostic(:w, OnGrid(Cells), u"kg"),
            Diagnostic(:k, OnGrid(Edges), u"J/s/m^3"),
            Prognostic(:x, OnGrid(Cells), u"J"),
        )
        model = Tile(strat, grid, DummyInitializer{:x}())
        # check vars
        @test hasproperty(model.state.vars,:top)
        @test hasproperty(model.state.vars,:testground1)
        @test hasproperty(model.state.vars,:testground2)
        @test hasproperty(model.state.vars,:bottom)
        # check for variables
        @test hasproperty(model.state.uproto, :x)
        @test hasproperty(model.state.uproto, :x)
        # clean-up method definitions (necessary for re-running test set)
        Base.delete_method(@which CryoGrid.variables(TestGroundLayer(TestGroundProcess()),TestGroundProcess()))
    end
    @testset "Helper functions" begin
        grid = Grid(Vector(0.0:10.0:1000.0)u"m")
        strat = Stratigraphy(
            -1.0u"m" => Top(TestBoundary()), (
                0.0u"m" => Named(:testgroundlayer1, TestGroundLayer(TestGroundProcess())),
                100.0u"m" => Named(:testgroundlayer2, TestGroundLayer(TestGroundProcess())),
            ),
            1000.0u"m" => Bottom(TestBoundary())
        )
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Prognostic(:x, OnGrid(Cells), u"J"),
            Diagnostic(:k, OnGrid(Edges), u"J/s/m^3"),
            Diagnostic(:w, OnGrid(Cells), u"kg"),
        )
        model = Tile(strat, grid, DummyInitializer{:x}())
        @test length(getvar(:x, model, model.state.uproto)) == length(cells(grid))
        @test length(getvar(:k, model, model.state.uproto)) == length(grid)
        # clean-up method definitions (necessary for re-running test set)
        Base.delete_method(@which CryoGrid.variables(TestGroundLayer(TestGroundProcess()),TestGroundProcess()))
    end
end
