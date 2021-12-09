module CryoGrid

global CRYOGRID_DEBUG = haskey(ENV,"CG_DEBUG") && ENV["CG_DEBUG"] == "true"
function debug(debug::Bool)
    @warn "Enabling/disabling debug mode after the module has been loaded will not update existing types and may cause errors."
    global CRYOGRID_DEBUG = debug
end

using Base: @propagate_inbounds
using Reexport

# Common types and methods
export Layer, SubSurface, Top, Bottom
export Process, SubSurfaceProcess, BoundaryProcess, CompoundProcess, Coupled
export BoundaryStyle, Dirichlet, Neumann
export Callback, CallbackStyle
include("types.jl")
export variables, initialcondition!, diagnosticstep!, prognosticstep!, interact!
export boundaryflux, boundaryvalue, criterion, affect!, observe
include("methods.jl")

# Submodules
include("Utils/Utils.jl")
include("Numerics/Numerics.jl")
include("Layers/Layers.jl")
include("Physics/Physics.jl")
include("Land/Land.jl")
include("IO/InputOutput.jl")
include("Diagnostics/Diagnostics.jl")
include("Drivers/Drivers.jl")

using .Utils
# Re-exported submodules
@reexport using .Numerics
@reexport using .Utils: convert_tspan
@reexport using .Layers
@reexport using .Physics
@reexport using .Land
@reexport using .Drivers
@reexport using .InputOutput
@reexport using .Diagnostics

# Re-exported packages
@reexport using Dates: Dates, DateTime
@reexport using DimensionalData: DimArray, Z, Dim, dims
@reexport using IntervalSets
@reexport using Unitful

# Include Presets submodule last to allow dependence on other submodules.
include("Presets/Presets.jl")

end # module
