using CryoGrid
using Test

include("types.jl")

@testset "Stratigraphy" begin
    @testset "3-layer" begin
        bounds = (-1.0u"m",0.0u"m",10.0u"m")
        strat = Stratigraphy(
            bounds[1] => Top(TestBoundary()),
            bounds[2] => Ground(:testground, TestGroundLayer(), TestGroundProcess()),
            bounds[3] => Bottom(TestBoundary())
        )
        # Check boundaries
        @test all([bounds[i] == b for (i,b) in enumerate(bounds)])
        # Check iteration and node types
        for (i,node) in enumerate(strat)
            if i == 1
                @test CryoGrid.nameof(node) == :top
                @test typeof(node.layer) <: Top
                @test typeof(node.process) <: Processes{Tuple{TestBoundary}}
            elseif i == 2
                @test CryoGrid.nameof(node) == :testground
                @test typeof(node.layer) <: TestGroundLayer
                @test typeof(node.process) <: Processes{Tuple{TestGroundProcess}}
            elseif i == 3
                @test CryoGrid.nameof(node) == :bottom
                @test typeof(node.layer) <: Bottom
                @test typeof(node.process) <: Processes{Tuple{TestBoundary}}
            end
        end
    end
    @testset "4-layer" begin
        bounds = (-1.0u"m",0.0u"m",10.0u"m",100.0u"m")
        strat = Stratigraphy(
            bounds[1] => Top(TestBoundary()),(
                bounds[2] => Ground(:testground1, TestGroundLayer(), TestGroundProcess()),
                bounds[3] => Ground(:testground2, TestGroundLayer(), TestGroundProcess())
            ),
            bounds[4] => Bottom(TestBoundary())
        )
        # Check boundaries
        @test all([bounds[i] == b for (i,b) in enumerate(bounds)])
        # Check iteration and node types
        for (i,node) in enumerate(strat)
            if i == 1
                @test CryoGrid.nameof(node) == :top
                @test typeof(node.layer) <: Top
                @test typeof(node.process) <: Processes{Tuple{TestBoundary}}
            elseif i == 2
                @test CryoGrid.nameof(node) == :testground1
                @test typeof(node.layer) <: TestGroundLayer
                @test typeof(node.process) <: Processes{Tuple{TestGroundProcess}}
            elseif i == 3
                @test CryoGrid.nameof(node) == :testground2
                @test typeof(node.layer) <: TestGroundLayer
                @test typeof(node.process) <: Processes{Tuple{TestGroundProcess}}
            elseif i == 4
                @test CryoGrid.nameof(node) == :bottom
                @test typeof(node.layer) <: Bottom
                @test typeof(node.process) <: Processes{Tuple{TestBoundary}}
            end
        end
    end
    bounds = (-1.0u"m",0.0u"m",10.0u"m",100.0u"m")
    # Check that duplicated layer names throws an error
    @test_throws AssertionError Stratigraphy(
        bounds[1] => Top(TestBoundary()), (
            bounds[2] => Ground(:duplicatedname, TestGroundLayer(), TestGroundProcess()),
            bounds[3] => Ground(:duplicatedname, TestGroundLayer(), TestGroundProcess()),
        ),
        bounds[4] => Bottom(TestBoundary())
    )
    # Check that mis-specified stratigraphies throw method errors
    @test_throws MethodError Stratigraphy(
        bounds[1] => Bottom(TestBoundary()),
        bounds[2] => Ground(:testground, TestGroundLayer(), TestGroundProcess()),
        bounds[3] => Top(TestBoundary())
    )
    @test_throws MethodError Stratigraphy(
        bounds[1] => Top(TestBoundary()),
        bounds[2] => Ground(:testground, TestGroundLayer(), TestGroundProcess()),
    )
    @test_throws MethodError Stratigraphy(
        bounds[1] => Ground(:testground, TestGroundLayer(), TestGroundProcess()),
        bounds[2] => Bottom(TestBoundary())
    )
end
