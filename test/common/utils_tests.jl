using CryoGrid.Utils: ffill!
using Test

@testset "Utils" begin
    @testset "ffill!" begin
        X = [missing, 1.0, 1.0, missing, missing, missing, 2.0, missing, 3.0, missing]
        ffill!(X)
        @test ismissing(X[1])
        @test X[2:end] == [1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0]
    end
end
