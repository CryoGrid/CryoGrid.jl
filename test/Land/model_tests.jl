using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils
using Test
using InteractiveUtils: @which

include("../types.jl")

@testset "3-layer" begin
    grid = Grid(Vector(0.0:10.0:1000.0)u"m")
    strat = Stratigraphy(
        -1.0u"m" => top(TestBoundary()),
        0.0u"m" => subsurface(:testground, TestGroundLayer(), TestGroundProcess()),
        1000.0u"m" => bottom(TestBoundary())
    )
    function checkfields(model)
        # for each non-error test case, we need to check that all layers are present in vars
        @test hasproperty(model.state.vars,:top)
        @test hasproperty(model.state.vars,:testground)
        @test hasproperty(model.state.vars,:bottom)
    end
    try
        # case: no variables defined
        @test_throws AssertionError model = LandModel(strat,grid)
        # case no prognostic variables defined
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
        )
        @test_throws AssertionError model = LandModel(strat,grid)
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Prognostic(:x, Float"J", OnGrid(Cells)),
            Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
        )
        model = LandModel(strat,grid)
        checkfields(model)
        @test hasproperty(model.state.uproto,:x)
        @test hasproperty(model.state.griddiag,:k)
        # warn when variable declared as prognostic and diagnostic
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Prognostic(:x, Float"J", OnGrid(Cells)),
            Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
            Diagnostic(:w, Float"kg", OnGrid(Cells)),
            Prognostic(:w, Float"kg", OnGrid(Cells)),
        )
        @test_logs (:warn,r".*declared as both prognostic/algebraic and diagnostic.*") model = LandModel(strat,grid)
        checkfields(model)
        @test hasproperty(model.state.uproto,:x)
        @test hasproperty(model.state.uproto,:w)
        # allow non-conflicting duplicates
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Prognostic(:x, Float"J", OnGrid(Cells)),
            Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
            Diagnostic(:w, Float"kg", OnGrid(Cells)),
            Diagnostic(:w, Float"kg", OnGrid(Cells)),
        )
        model = LandModel(strat,grid)
        checkfields(model)
        @test hasproperty(model.state.griddiag,:w)
        # error when conflicting variables declared
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Prognostic(:x, Float"J", OnGrid(Cells)),
            Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
            Diagnostic(:w, Float"kg/m", OnGrid(Cells)),
            Diagnostic(:w, Float"kg/m", OnGrid(Edges)),
        )
        @test_throws AssertionError model = LandModel(strat,grid)
        # test scalar variables and parameters
        CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
            Prognostic(:x, Float"J", OnGrid(Cells)),
            Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
            Diagnostic(:a, Float64, Scalar),
        )
        model = LandModel(strat,grid)
        checkfields(model)
        @test hasproperty(model.state.diag.testground, :a)
        model.state.diag.testground.a.cache.du[1] = 2.0
        @test model.state.diag.testground.a.cache.du[1] == 2.0
        state = getstate(:testground, model, model.state.uproto, model.state.uproto, 0.0)
    finally
        # clean-up method definitions (necessary for re-running test set)
        Base.delete_method(@which CryoGrid.variables(TestGroundLayer(),TestGroundProcess()))
    end
end
@testset "4-layer" begin
    grid = Grid(Vector(0.0:10.0:1000.0)u"m")
    strat = Stratigraphy(
        -1.0u"m" => top(TestBoundary()), (
            0.0u"m" => subsurface(:testground1, TestGroundLayer(), TestGroundProcess()),
            100.0u"m" => subsurface(:testground2, TestGroundLayer(), TestGroundProcess()),
        ),
        1000.0u"m" => bottom(TestBoundary())
    )
    CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
        Diagnostic(:w,Float"kg",OnGrid(Cells)),
        Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
        Prognostic(:x, Float"J", OnGrid(Cells)),
    )
    model = LandModel(strat,grid)
    # check vars
    @test hasproperty(model.state.vars,:top)
    @test hasproperty(model.state.vars,:testground1)
    @test hasproperty(model.state.vars,:testground2)
    @test hasproperty(model.state.vars,:bottom)
    # check for variables
    @test hasproperty(model.state.uproto, :x)
    @test hasproperty(model.state.uproto, :x)
    # clean-up method definitions (necessary for re-running test set)
    Base.delete_method(@which CryoGrid.variables(TestGroundLayer(),TestGroundProcess()))
end
@testset "Helper functions" begin
    grid = Grid(Vector(0.0:10.0:1000.0)u"m")
    strat = Stratigraphy(
        -1.0u"m" => top(TestBoundary()), (
            0.0u"m" => subsurface(:testground1, TestGroundLayer(), TestGroundProcess()),
            100.0u"m" => subsurface(:testground2, TestGroundLayer(), TestGroundProcess()),
        ),
        1000.0u"m" => bottom(TestBoundary())
    )
    CryoGrid.variables(::TestGroundLayer, ::TestGroundProcess) = (
        Prognostic(:x, Float"J", OnGrid(Cells)),
        Diagnostic(:k, Float"J/s/m^3", OnGrid(Edges)),
        Diagnostic(:w,Float"kg",OnGrid(Cells)),
    )
    model = LandModel(strat,grid)
    @test length(getvar(:x, model, model.state.uproto)) == length(cells(grid))
    @test length(getvar(:k, model, model.state.uproto)) == length(grid)
    # clean-up method definitions (necessary for re-running test set)
    Base.delete_method(@which CryoGrid.variables(TestGroundLayer(),TestGroundProcess()))
end
