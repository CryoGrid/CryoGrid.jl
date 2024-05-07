using CryoGrid
using Test

@testset "Stratigraphy" begin
    @testset "3-layer" begin
        bounds = (-1.0u"m",0.0u"m",10.0u"m")
        strat = Stratigraphy(
            bounds[1] => Top(TestBoundary()),
            bounds[2] => TestGroundLayer(TestGroundProcess()),
            bounds[3] => Bottom(TestBoundary())
        )
        # Check boundaries
        @test all([bounds[i] == b for (i,b) in enumerate(bounds)])
        # Check iteration and layer types
        for (i,named_layer) in enumerate(namedlayers(strat))
            if i == 1
                @test nameof(named_layer) == :top
                @test typeof(named_layer.val) <: Top
                @test typeof(named_layer.val.proc) <: TestBoundary
            elseif i == 2
                @test nameof(named_layer) == :testgroundlayer
                @test typeof(named_layer.val) <: TestGroundLayer
                @test typeof(named_layer.val.proc) <: TestGroundProcess
            elseif i == 3
                @test nameof(named_layer) == :bottom
                @test typeof(named_layer.val) <: Bottom
                @test typeof(named_layer.val.proc) <: TestBoundary
            end
        end
    end
    @testset "4-layer" begin
        bounds = (-1.0u"m",0.0u"m",10.0u"m",100.0u"m")
        strat = Stratigraphy(
            bounds[1] => Top(TestBoundary()),(
                bounds[2] => TestGroundLayer(TestGroundProcess()),
                bounds[3] => TestGroundLayer(TestGroundProcess()),
            ),
            bounds[4] => Bottom(TestBoundary())
        )
        # Check boundaries
        @test all([bounds[i] == b for (i,b) in enumerate(bounds)])
        # Check iteration and layer types
        for (i,named_layer) in enumerate(namedlayers(strat))
            if i == 1
                @test nameof(named_layer) == :top
                @test typeof(named_layer.val) <: Top
                @test typeof(named_layer.val.proc) <: TestBoundary
            elseif i == 2
                @test nameof(named_layer) == :testgroundlayer1
                @test typeof(named_layer.val) <: TestGroundLayer
                @test typeof(named_layer.val.proc) <: TestGroundProcess
            elseif i == 3
                @test nameof(named_layer) == :testgroundlayer2
                @test typeof(named_layer.val) <: TestGroundLayer
                @test typeof(named_layer.val.proc) <: TestGroundProcess
            elseif i == 4
                @test nameof(named_layer) == :bottom
                @test typeof(named_layer.val) <: Bottom
                @test typeof(named_layer.val.proc) <: TestBoundary
            end
        end
    end
    @testset "4-layer w/ coupling" begin
        bounds = (-1.0u"m",0.0u"m",10.0u"m",100.0u"m")
        strat = Stratigraphy(
            bounds[1] => Top(Coupled(TestBoundary(), TestBoundary())),(
                bounds[2] => TestGroundLayer(Coupled(TestGroundProcess(), TestGroundProcess())),
                bounds[3] => TestGroundLayer(TestGroundProcess()),
            ),
            bounds[4] => Bottom(Coupled(TestBoundary(), TestBoundary()))
        )
        # Check boundaries
        @test all([bounds[i] == b for (i,b) in enumerate(bounds)])
        # Check iteration and component types
        for (i,named_layer) in enumerate(namedlayers(strat))
            if i == 1
                @test nameof(named_layer) == :top
                @test typeof(named_layer.val) <: Top
                @test typeof(named_layer.val.proc) <: CoupledProcesses{Tuple{TestBoundary,TestBoundary}}
            elseif i == 2
                @test nameof(named_layer) == :testgroundlayer1
                @test typeof(named_layer.val) <: TestGroundLayer
                @test typeof(named_layer.val.proc) <: CoupledProcesses{Tuple{TestGroundProcess,TestGroundProcess}}
            elseif i == 3
                @test nameof(named_layer) == :testgroundlayer2
                @test typeof(named_layer.val) <: TestGroundLayer
                @test typeof(named_layer.val.proc) <: TestGroundProcess
            elseif i == 4
                @test nameof(named_layer) == :bottom
                @test typeof(named_layer.val) <: Bottom
                @test typeof(named_layer.val.proc) <: CoupledProcesses{Tuple{TestBoundary,TestBoundary}}
            end
        end
    end
    bounds = (-1.0u"m",0.0u"m",10.0u"m",100.0u"m")
    # Check that layer names are automatically deduplicated
    strat = Stratigraphy(
        bounds[1] => Top(TestBoundary()), (
            bounds[2] => Named(:duplicatedname, TestGroundLayer(TestGroundProcess())),
            bounds[3] => Named(:duplicatedname, TestGroundLayer(TestGroundProcess())),
        ),
        bounds[4] => Bottom(TestBoundary())
    )
    @test hasproperty(strat, :duplicatedname1)
    @test hasproperty(strat, :duplicatedname2)
    # Check that mis-specified stratigraphies throw method errors
    @test_throws MethodError Stratigraphy(
        bounds[1] => Bottom(TestBoundary()),
        bounds[2] => TestGroundLayer(TestGroundProcess()),
        bounds[3] => Top(TestBoundary())
    )
    @test_throws MethodError Stratigraphy(
        bounds[1] => Top(TestBoundary()),
        bounds[2] => TestGroundLayer(TestGroundProcess()),
    )
    @test_throws MethodError Stratigraphy(
        bounds[1] => TestGroundLayer(TestGroundProcess()),
        bounds[2] => Bottom(TestBoundary())
    )
    # Check that out-of-order bounds throw an error
    @test_throws AssertionError Stratigraphy(
        bounds[2] => Top(TestBoundary()), (
            bounds[3] => TestGroundLayer(TestGroundProcess()),
            bounds[1] => TestGroundLayer(TestGroundProcess()),
        ),
        bounds[4] => Bottom(TestBoundary())
    )
    # Check that duplicated boundaries are allowed
    bounds = (0.0u"m",0.0u"m",10.0u"m")
    @test_nowarn strat = Stratigraphy(
        bounds[1] => Top(TestBoundary()),
        bounds[2] => TestGroundLayer(TestGroundProcess()),
        bounds[3] => Bottom(TestBoundary())
    )
end
