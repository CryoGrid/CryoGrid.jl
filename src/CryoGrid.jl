module CryoGrid

global CRYOGRID_DEBUG = haskey(ENV,"CG_DEBUG") && ENV["CG_DEBUG"] == "true"
function debug(debug::Bool)
    @warn "Enabling/disabling debug mode after the module has been loaded will not update existing types and may cause errors."
    global CRYOGRID_DEBUG = debug
end

using Base: @propagate_inbounds
using Reexport

export Layer, SubSurface, Top, Bottom, Boundary
export Process, SubSurfaceProcess, BoundaryProcess, CompositeProcess, Coupled
export BoundaryStyle, Dirichlet, Neumann
export AbstractParameterization, Parameterization
export variables, initialcondition!, diagnosticstep!, prognosticstep!, interact!, observe

# Common types and methods
include("Common/types.jl")
include("Common/methods.jl")

# Submodules
include("Common/Utils/Utils.jl")
include("Common/Numerics/Numerics.jl")
include("Common/Forcings/Forcings.jl")
include("IO/InputOutput.jl")
include("Layers/Layers.jl")
include("Processes/Processes.jl")
include("Setup/Setup.jl")
include("Diagnostics/Diagnostics.jl")
include("Callbacks/Callbacks.jl")

# Re-exported submodules
@reexport using .Utils
@reexport using .Numerics
@reexport using .Forcings
@reexport using .InputOutput: loadforcings, InputSpec, JsonSpec
@reexport using .Layers
@reexport using .Processes
@reexport using .Setup
@reexport using .Diagnostics
@reexport using .Callbacks

# Re-exported packages
@reexport using Dates: Dates, DateTime
@reexport using DimensionalData: DimArray, Z, Dim, dims
@reexport using IntervalSets
@reexport using Unitful

# Import parameters function into top-level scope;
# We do not export it in order to avoid naming conflicts with other packages.
import .Setup: parameters

# Include Models submodule last to allow dependence on other submodules.
include("Models/Models.jl")

end # module
