global CRYOGRID_DEBUG = haskey(ENV,"CG_DEBUG") && ENV["CG_DEBUG"] == "true"

include("utils.jl")
include("types.jl")
include("math.jl")
include("grid.jl")
include("forcing.jl")
include("variables.jl")
include("stratigraphy.jl")
include("setup.jl")
include("problem.jl")
include("output.jl")
