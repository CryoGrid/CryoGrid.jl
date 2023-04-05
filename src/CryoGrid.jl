module CryoGrid

global CRYOGRID_DEBUG = haskey(ENV,"CG_DEBUG") && ENV["CG_DEBUG"] == "true"
function debug(debug::Bool)
    global CRYOGRID_DEBUG = debug
    # disable loop vectorization in debug mode
    Numerics.turbo(!debug)
    CRYOGRID_DEBUG && @warn "Debug mode enabled! Some performance features such as loop vectorization are now turned off by default."
    return CRYOGRID_DEBUG
end

using Base: @propagate_inbounds
using ComponentArrays
using ConstructionBase
using Dates
using IfElse
using LinearAlgebra
using ModelParameters
using Reexport

import Flatten
import Interpolations

# Re-exported packages
@reexport using Dates: Dates, Date, DateTime
@reexport using DimensionalData
@reexport using IfElse
@reexport using IntervalSets
@reexport using ModelParameters
@reexport using Setfield: @set, @set!
@reexport using Unitful
@reexport using UnPack

export Interpolations

# Common types and methods
export Layer, SubSurface, Top, Bottom
export Process, SubSurfaceProcess, BoundaryProcess, CoupledProcesses
export Coupled, Coupled2, Coupled3, Coupled4
export DiscreteEvent, ContinuousEvent, GridContinuousEvent
include("types.jl")
export BoundaryCondition, hasfixedvolume
include("traits.jl")
export variables, basevariables, processes, initialcondition!, diagnosticstep!, prognosticstep!, interact!, timestep
export boundaryflux, boundaryvalue, criterion, criterion!, trigger!
include("methods.jl")

# Submodules
include("Utils/Utils.jl")
using .Utils
export convert_t, convert_tspan, pstrip, @pstrip, @sym_str
include("Numerics/Numerics.jl")
using .Numerics
export DiscretizationStrategy, AutoGrid, PresetGrid, Grid, cells, edges, subgridinds, Δ, volume, area
export initializer, getvar
include("IO/InputOutput.jl")
@reexport using .InputOutput
include("Physics/Physics.jl")
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

# include dependent submodules
include("Solvers/Solvers.jl")
include("Presets/Presets.jl")

end # module
