using CryoGrid
using CryoGrid.Numerics
using Test, BenchmarkTools

include("../testutils.jl")

@testset "Grid" begin
    gridvals = [1:0.1:2...,2.5:0.5:3...,4:1.0:10...]u"m"
    grid = Grid(gridvals)
    @testset "Functions" begin
        @test all([grid[i] == gridvals[i] for (i,val) in enumerate(gridvals)])
        cgrid = @inferred cells(grid)
        @test all([cgrid[i] ≈ (grid[i]+grid[i+1])/2 for i in 1:length(cgrid)])
        egrid = @inferred edges(cgrid)
        @test all([egrid[i] ≈ grid[i] for i in 1:length(grid)])
        Δgrid = @inferred Δ(grid)
        @test all([Δgrid[i-1] ≈ grid[i] - grid[i-1] for i in 2:length(grid)])
        Δcgrid = @inferred Δ(cgrid)
        @test all([Δcgrid[i-1] ≈ cgrid[i] - cgrid[i-1] for i in 2:length(cgrid)])
        @test all(grid[1.0u"m"..10.0u"m"] .≈ grid)
        @test all(grid[1.0u"m"..2.0u"m"] .≈ grid[1:11])
        @test all(grid[Interval{:open,:open}(1.0u"m",10.0u"m")] .≈ grid[2:end-1])
        @test all(grid[Interval{:closed,:open}(1.0u"m",10.0u"m")] .≈ grid[1:end-1])
        @test all(grid[Interval{:open,:closed}(1.0u"m",10.0u"m")] .≈ grid[2:end])
    end
    @testset "Allocations" begin
        benchres = @benchmark Δ($grid)
        @test benchres.allocs == 0
    end
end
