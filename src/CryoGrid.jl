module CryoGrid

global CRYOGRID_DEBUG = haskey(ENV,"CG_DEBUG") && ENV["CG_DEBUG"] == "true"
function debug(debug::Bool)
    @warn "Enabling/disabling debug mode after the module has been loaded will not update existing types and may cause errors."
    global CRYOGRID_DEBUG = debug
end

using Base: @propagate_inbounds
using ComponentArrays
using Dates
using LinearAlgebra
using Reexport

# Common types and methods
export Layer, SubSurface, Top, Bottom
export Process, SubSurfaceProcess, BoundaryProcess, CoupledProcesses
export DiscreteEvent, ContinuousEvent, GridContinuousEvent
export Coupled, Coupled2, Coupled3, Coupled4
include("types.jl")
export variables, basevariables, processes, initialcondition!, diagnosticstep!, prognosticstep!, interact!, timestep
export boundaryflux, boundaryvalue, criterion, criterion!, trigger!, observe
include("methods.jl")

# Submodules
include("Utils/Utils.jl")
using .Utils
export convert_t, convert_tspan, pstrip, @pstrip, @sym_str
include("Numerics/Numerics.jl")
using .Numerics
export Grid, cells, edges, subgridinds, Δ, volume, area, initializer, getvar
include("IO/InputOutput.jl")
@reexport using .InputOutput
include("Physics/Physics.jl")
using .Physics
include("Strat/Strat.jl")
@reexport using .Strat
parameters = Strat.parameters
include("Diagnostics/Diagnostics.jl")
@reexport using .Diagnostics

# Coupling
include("coupling.jl")

# Problem interface
export CryoGridProblem
include("problem.jl")

# Re-exported packages
@reexport using Dates: Dates, DateTime
@reexport using DimensionalData
@reexport using IntervalSets
@reexport using Unitful

# include dependent submodules
include("Solvers/Solvers.jl")
include("Presets/Presets.jl")

end # module
