using CryoGrid
using Dates
using Test

struct TestCallable end
(::TestCallable)(x,y,z) = nothing

@testset "Utils" begin
    @testset "argnames" begin
        f(x) = nothing
        g(x,y) = nothing
        h(x,y,args...) = nothing
        @test Utils.argnames(f) == [:x,]
        @test Utils.argnames(g) == [:x,:y]
        @test Utils.argnames(h) == [:x,:y,:args]
        @test Utils.argnames(TestCallable()) == [:x,:y,:z]
    end
    @testset "convert time" begin
        epoch_date = Dates.epochms2datetime(0.0)
        # solver time is seconds since epoch, so we check that the conversion function does it correctly
        @test convert_t(epoch_date + Hour(1)) ≈ Dates.value(convert(Second, Hour(1)))
        @test convert_t(convert_t(epoch_date + Hour(1))) == epoch_date + Hour(1)
        @test all(convert_tspan((epoch_date, epoch_date + Hour(1))) .≈ (0.0, Dates.value(convert(Second, Hour(1)))))
        @test all(convert_tspan(convert_tspan((epoch_date, epoch_date + Hour(1)))) .== (epoch_date, epoch_date + Hour(1)))
    end
    @testset "ffill!" begin
        X = [missing, 1.0, 1.0, missing, missing, missing, 2.0, missing, 3.0, missing]
        Utils.ffill!(X)
        @test ismissing(X[1])
        @test X[2:end] == [1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0]
    end
    @testset "getscalar" begin
        @test getscalar([1.0]) == 1.0
        @test getscalar(1.0) == 1.0
        @test getscalar([1.0,2.0], 2) == 2.0
        @test_throws AssertionError getscalar([1.0,2.0])
    end
    @testset "tuplejoin" begin
        @test tuplejoin(()) == ()
        @test tuplejoin((),(1,)) == (1,)
        @test tuplejoin(1,2,3) == (1,2,3)
        @test tuplejoin((1,),(2,),(3,)) == (1,2,3)
        @test tuplejoin((1,2,3),(4,5)) == (1,2,3,4,5)
        # should not flatten nested tuples
        @test tuplejoin((1,2,(3,)),(4,)) == (1,2,(3,),4)
    end
end
