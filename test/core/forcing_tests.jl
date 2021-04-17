using CryoGrid
using TimeSeries, Dates
using Test, BenchmarkTools

@testset "TimeSeriesForcing" begin
    ts = DateTime(2020,1,1):Day(1):DateTime(2020,3,1)
    ys = randn(length(ts))
    forcing = CryoGrid.TimeSeriesForcing(ys, ts, :test)
    @test all(timestamp(forcing.tarray) .== ts)
    @test forcing(Dates.datetime2epochms(ts[1])/1000.0) ≈ ys[1]
    @test forcing(Dates.datetime2epochms(ts[end])/1000.0) ≈ ys[end]
    @test_throws BoundsError forcing(-1.0)
    @test_throws BoundsError forcing(Dates.datetime2epochms(ts[end])/1000.0+1.0)
    t1,y1 = ts[1], ys[1]
    t2,y2 = ts[2], ys[2]
    @test forcing((Dates.datetime2epochms(t1) + Dates.datetime2epochms(t2))/2000.0) ≈ (y1+y2)/2
    t = Dates.datetime2epochms(t1)/1000.0
    benchres = @benchmark $forcing($t)
    @test benchres.allocs == 0
    out = zeros(100)
    queries = t .+ (1:100);
    benchres = @benchmark $out .= $forcing.($queries)
    @test benchres.allocs == 0
end
