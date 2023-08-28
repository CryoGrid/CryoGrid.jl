@testset "Solvers" begin
    include("../test_problems.jl")
    include("cgeuler_tests.jl")
    include("cglite_implicit_tests.jl")
end
