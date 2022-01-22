using CryoGrid

@testset "Parameters" begin
    @testset "LinearTrend" begin
       trend = LinearTrend(slope=0.0, intercept=1.0)
       @test CryoGrid.Strat.transform((t=1.0,), trend) ≈ 1.0
       trend = LinearTrend(slope=0.1, intercept=1.0)
       @test CryoGrid.Strat.transform((t=1.0,), trend) ≈ 1.1
       trend = LinearTrend(slope=0.1, intercept=0.5, minval=0.0, maxval=1.0)
       @test CryoGrid.Strat.transform((t=10.0,), trend) ≈ 1.0
       trend = LinearTrend(slope=-0.1, intercept=0.5, minval=0.0, maxval=1.0)
       @test CryoGrid.Strat.transform((t=10.0,), trend) ≈ 0.0
       trend = LinearTrend(slope=0.1, intercept=0.5, tstart=0.0, tstop=1.0)
       @test CryoGrid.Strat.transform((t=2.0,), trend) ≈ 0.6
    end
    @testset "PiecewiseLinear" begin
        pc = PiecewiseLinear(0.0, (1.0,1.0); tstart=0.0, tstop=1.0)
        @test CryoGrid.Strat.transform((t=-0.1,), pc) ≈ 0.0
        @test CryoGrid.Strat.transform((t=0.0,), pc) ≈ 0.0
        @test CryoGrid.Strat.transform((t=0.5,), pc) ≈ 0.5
        @test CryoGrid.Strat.transform((t=1.0,), pc) ≈ 1.0
        @test CryoGrid.Strat.transform((t=1.1,), pc) ≈ 1.0
        pc = PiecewiseLinear(1.0, (0.4,0.5), (0.6,0.0); tstart=0.0, tstop=1.0)
        @test CryoGrid.Strat.transform((t=-0.1,), pc) ≈ 1.0
        @test CryoGrid.Strat.transform((t=0.0,), pc) ≈ 1.0
        @test CryoGrid.Strat.transform((t=0.7,), pc) ≈ 0.25
        @test CryoGrid.Strat.transform((t=1.0,), pc) ≈ 0.0
    end
end
