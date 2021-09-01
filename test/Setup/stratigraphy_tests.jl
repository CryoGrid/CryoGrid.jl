using CryoGrid
using Test

include("../types.jl")

@testset "Stratigraphy" begin
    @testset "3-layer" begin
        bounds = (-1.0u"m",0.0u"m",10.0u"m")
        strat = Stratigraphy(
            bounds[1] => top(TestBoundary()),
            bounds[2] => subsurface(:testground, TestGroundLayer(), TestGroundProcess()),
            bounds[3] => bottom(TestBoundary())
        )
        # Check boundaries
        @test all([bounds[i] == b for (i,b) in enumerate(bounds)])
        # Check iteration and node types
        for (i,node) in enumerate(strat)
            if i == 1
                @test CryoGrid.nodename(node) == :top
                @test typeof(node.layer) <: Top
                @test typeof(node.process) <: System{Tuple{TestBoundary}}
            elseif i == 2
                @test CryoGrid.nodename(node) == :testground
                @test typeof(node.layer) <: TestGroundLayer
                @test typeof(node.process) <: System{Tuple{TestGroundProcess}}
            elseif i == 3
                @test CryoGrid.nodename(node) == :bottom
                @test typeof(node.layer) <: Bottom
                @test typeof(node.process) <: System{Tuple{TestBoundary}}
            end
        end
    end
    @testset "4-layer" begin
        bounds = (-1.0u"m",0.0u"m",10.0u"m",100.0u"m")
        strat = Stratigraphy(
            bounds[1] => top(TestBoundary()),(
                bounds[2] => subsurface(:testground1, TestGroundLayer(), TestGroundProcess()),
                bounds[3] => subsurface(:testground2, TestGroundLayer(), TestGroundProcess())
            ),
            bounds[4] => bottom(TestBoundary())
        )
        # Check boundaries
        @test all([bounds[i] == b for (i,b) in enumerate(bounds)])
        # Check iteration and node types
        for (i,node) in enumerate(strat)
            if i == 1
                @test CryoGrid.nodename(node) == :top
                @test typeof(node.layer) <: Top
                @test typeof(node.process) <: System{Tuple{TestBoundary}}
            elseif i == 2
                @test CryoGrid.nodename(node) == :testground1
                @test typeof(node.layer) <: TestGroundLayer
                @test typeof(node.process) <: System{Tuple{TestGroundProcess}}
            elseif i == 3
                @test CryoGrid.nodename(node) == :testground2
                @test typeof(node.layer) <: TestGroundLayer
                @test typeof(node.process) <: System{Tuple{TestGroundProcess}}
            elseif i == 4
                @test CryoGrid.nodename(node) == :bottom
                @test typeof(node.layer) <: Bottom
                @test typeof(node.process) <: System{Tuple{TestBoundary}}
            end
        end
    end
    bounds = (-1.0u"m",0.0u"m",10.0u"m",100.0u"m")
    # Check that duplicated layer names throws an error
    @test_throws AssertionError Stratigraphy(
        bounds[1] => top(TestBoundary()), (
            bounds[2] => subsurface(:duplicatedname, TestGroundLayer(), TestGroundProcess()),
            bounds[3] => subsurface(:duplicatedname, TestGroundLayer(), TestGroundProcess()),
        ),
        bounds[4] => bottom(TestBoundary())
    )
    # Check that mis-specified stratigraphies throw method errors
    @test_throws MethodError Stratigraphy(
        bounds[1] => bottom(TestBoundary()),
        bounds[2] => subsurface(:testground, TestGroundLayer(), TestGroundProcess()),
        bounds[3] => top(TestBoundary())
    )
    @test_throws MethodError Stratigraphy(
        bounds[1] => top(TestBoundary()),
        bounds[2] => subsurface(:testground, TestGroundLayer(), TestGroundProcess()),
    )
    @test_throws MethodError Stratigraphy(
        bounds[1] => subsurface(:testground, TestGroundLayer(), TestGroundProcess()),
        bounds[2] => bottom(TestBoundary())
    )
    # Check that out-of-order bounds throw an error
    @test_throws AssertionError Stratigraphy(
        bounds[2] => top(TestBoundary()), (
            bounds[3] => subsurface(:duplicatedname, TestGroundLayer(), TestGroundProcess()),
            bounds[1] => subsurface(:duplicatedname, TestGroundLayer(), TestGroundProcess()),
        ),
        bounds[4] => bottom(TestBoundary())
    )
    # Check that duplicated boundaries are allowed
    bounds = (0.0u"m",0.0u"m",10.0u"m")
    @test_nowarn strat = Stratigraphy(
        bounds[1] => top(TestBoundary()),
        bounds[2] => subsurface(:testground, TestGroundLayer(), TestGroundProcess()),
        bounds[3] => bottom(TestBoundary())
    )
end
