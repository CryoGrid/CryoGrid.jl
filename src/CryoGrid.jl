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
export Process, SubSurfaceProcess, BoundaryProcess, CoupledProcesses, Coupled
include("types.jl")
export variables, initialcondition!, diagnosticstep!, prognosticstep!, interact!
export boundaryflux, boundaryvalue, criterion, affect!, observe
include("methods.jl")

# Submodules
include("Utils/Utils.jl")
include("Numerics/Numerics.jl")
include("IO/InputOutput.jl")
include("Physics/Physics.jl")
include("Strat/Strat.jl")
include("Diagnostics/Diagnostics.jl")

using .Numerics
export Grid, cells, edges, subgridinds, Δ, volume, area, initializer, getvar
using .Utils
export convert_t, convert_tspan, deparam, @sym_str
# Re-exported submodules
@reexport using .Physics
@reexport using .Strat
@reexport using .InputOutput
@reexport using .Diagnostics

# Re-exported packages
@reexport using Dates: Dates, DateTime
@reexport using DimensionalData: DimArray, Z, Dim, dims
@reexport using IntervalSets
@reexport using Unitful

# include dependent submodules
include("Drivers/Drivers.jl")
@reexport using .Drivers
include("Presets/Presets.jl")

end # module
