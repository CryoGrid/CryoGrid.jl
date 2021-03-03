@testset "Core" begin
    include("math_tests.jl")
    include("forcing_tests.jl")
    include("grid_tests.jl")
    include("stratigraphy_tests.jl")
    include("setup_tests.jl")
end
