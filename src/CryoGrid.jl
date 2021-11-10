module CryoGrid

global CRYOGRID_DEBUG = haskey(ENV,"CG_DEBUG") && ENV["CG_DEBUG"] == "true"
function debug(debug::Bool)
    @warn "Enabling/disabling debug mode after the module has been loaded will not update existing types and may cause errors."
    global CRYOGRID_DEBUG = debug
end

using Base: @propagate_inbounds
using Reexport

export Layer, SubSurface, Top, Bottom, Boundary
export Process, SubSurfaceProcess, BoundaryProcess, CompoundProcess, Coupled
export BoundaryStyle, Dirichlet, Neumann
export AbstractParameterization, Parameterization
export variables, initialcondition!, diagnosticstep!, prognosticstep!, interact!, observe

# Common types and methods
include("types.jl")
include("methods.jl")

# Submodules
include("Utils/Utils.jl")
include("Numerics/Numerics.jl")
include("IO/InputOutput.jl")
include("Layers/Layers.jl")
include("Physics/Physics.jl")
include("Land/Land.jl")
include("Drivers/Drivers.jl")
include("Diagnostics/Diagnostics.jl")
include("Callbacks/Callbacks.jl")

# Re-exported submodules
@reexport using .Utils
@reexport using .Numerics
@reexport using .InputOutput
@reexport using .Layers
@reexport using .Physics
@reexport using .Land
@reexport using .Drivers
@reexport using .Diagnostics
@reexport using .Callbacks

# Re-exported packages
@reexport using Dates: Dates, DateTime
@reexport using DimensionalData: DimArray, Z, Dim, dims
@reexport using IntervalSets
@reexport using Unitful

# Include Presets submodule last to allow dependence on other submodules.
include("Presets/Presets.jl")

end # module
