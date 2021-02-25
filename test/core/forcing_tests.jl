using CryoGrid
using TimeSeries, Dates
using Test

@testset "TimeSeriesForcing" begin
    ts = DateTime(2020,1,1):Day(1):DateTime(2020,3,1)
    ys = randn(length(ts))
    forcing = CryoGrid.TimeSeriesForcing(ys, ts, :test)
    @test forcing.t0 == ts[1]
    @test all(timestamp(forcing.tarray) .== ts)
    @test forcing(0.0) ≈ ys[1]
    @test forcing(Dates.value(ts[end]-ts[1])/1000.0) ≈ ys[end]
    @test_throws BoundsError forcing(-1.0)
    @test_throws BoundsError forcing(1.0e10)
    t1,y1 = ts[1], ys[1]
    t2,y2 = ts[2], ys[2]
    @test forcing(Dates.value((t2-t1)/2)/1000.0) ≈ (y1+y2)/2
end
