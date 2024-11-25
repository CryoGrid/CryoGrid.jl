using CryoGrid
using Dates
using Test, BenchmarkTools

@testset "Interpolated1D" begin
    ts = DateTime(2020,1,1):Day(1):DateTime(2020,3,1)
    ys = randn(length(ts))
    arr = DimStack((test=DimArray(ys, (Ti(ts),)),))
    interp = InputOutput.Interpolated1D(arr)
    f = interp.test
    @test f(Dates.datetime2epochms(ts[1])/1000.0) ≈ ys[1]
    @test f(Dates.datetime2epochms(ts[end])/1000.0) ≈ ys[end]
    @test_throws BoundsError f(-1.0)
    @test_throws BoundsError f(Dates.datetime2epochms(ts[end])/1000.0+1.0)
    t1,y1 = ts[1], ys[1]
    t2,y2 = ts[2], ys[2]
    @test f((Dates.datetime2epochms(t1) + Dates.datetime2epochms(t2))/2000.0) ≈ (y1+y2)/2
    t = Dates.datetime2epochms(t1)/1000.0
    benchres = @benchmark $f($t)
    @test benchres.allocs == 0
    out = zeros(100)
    queries = t .+ (1:100);
    benchres = @benchmark $out .= $f.($queries)
    @test benchres.allocs == 0
end