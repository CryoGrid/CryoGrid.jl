using CryoGrid
using Test, BenchmarkTools

@testset "Grid" begin
    gridvals = [1:0.1:2...,2.5:0.5:3...,4:1.0:10...]u"m"
    grid = Grid(gridvals)
    @test all([grid[i] == gridvals[i] for (i,val) in enumerate(gridvals)])
    cgrid = cells(grid)
    @test all([cgrid[i] ≈ (grid[i]+grid[i+1])/2 for i in 1:length(cgrid)])
    egrid = edges(cgrid)
    @test all([egrid[i] ≈ grid[i] for i in 1:length(grid)])
    Δgrid = Δ(grid)
    @test all([Δgrid[i-1] ≈ grid[i] - grid[i-1] for i in 2:length(grid)])
    Δcgrid = Δ(cgrid)
    @test all([Δcgrid[i-1] ≈ cgrid[i] - cgrid[i-1] for i in 2:length(cgrid)])
    @test all(grid[1.0u"m"..10.0u"m"] .≈ grid)
    @test all(grid[1.0u"m"..2.0u"m"] .≈ grid[1:11])
    @test all(grid[Interval{:open,:open}(1.0u"m",10.0u"m")] .≈ grid[2:end-1])
    @test all(grid[Interval{:closed,:open}(1.0u"m",10.0u"m")] .≈ grid[1:end-1])
    @test all(grid[Interval{:open,:closed}(1.0u"m",10.0u"m")] .≈ grid[2:end])
    @testset "allocations" begin
        benchres = @benchmark Δ($grid)
        @test benchres.allocs == 0
    end
    @testset "type stability" begin
        @inferred cells(grid)
        @inferred edges(cgrid)
        @inferred Δ(grid)
        # Note: constructor (and thus also subgridding) are *NOT* type stable
    end
end
