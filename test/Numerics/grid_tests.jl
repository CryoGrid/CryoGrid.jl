using CryoGrid
using CryoGrid.Numerics
using Test, BenchmarkTools

@testset "Grid" begin
    gridvals = [1:0.1:2...,2.5:0.5:3...,4:1.0:10...]u"m"
    grid = Grid(gridvals)
    @testset "Functions" begin
        @test length(grid) == length(gridvals)
        @test size(grid) == size(gridvals)
        @test length(cells(grid)) == length(gridvals) - 1
        @test firstindex(grid) == 1
        @test lastindex(grid) == length(grid)
        @test firstindex(grid[1..10]) == 1
        @test lastindex(grid[1..10]) == 10
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
        @test all(grid[1..10] .≈ grid[1:10])
        @test all(grid[5..10][2..4] .≈ grid[6:8])
        @test all(grid[1..10] .≈ grid[1:10])
        @test all(grid[Interval{:open,:open}(1.0u"m",10.0u"m")] .≈ grid[2:end-1])
        @test all(grid[Interval{:closed,:open}(1.0u"m",10.0u"m")] .≈ grid[1:end-1])
        @test all(grid[Interval{:open,:closed}(1.0u"m",10.0u"m")] .≈ grid[2:end])
        @test parent(grid[2..10]) == grid
        newgrid = similar(grid)
        updategrid!(newgrid, grid.*2.0u"m/m" .+ 1.0u"m")
        @test all(newgrid .≈ grid.*2.0u"m/m" .+ 1.0u"m")
        @test all(cells(newgrid) .≈ cells(grid).*2.0u"m/m" .+ 1.0u"m")
        @test all(Δ(newgrid) .≈ Δ(grid).*2.0u"m/m")
        @test all(Δ(cells(newgrid)) .≈ Δ(cells(grid)).*2.0u"m/m")
    end
    @testset "Allocations" begin
        benchres = @benchmark Δ($grid)
        @test benchres.allocs == 0
    end
end
