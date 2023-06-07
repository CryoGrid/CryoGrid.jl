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
@reexport using Dates
@reexport using DiffEqBase
@reexport using DiffEqCallbacks
@reexport using DimensionalData
@reexport using IfElse
@reexport using IntervalSets
@reexport using ModelParameters
@reexport using Setfield: @set, @set!
@reexport using SciMLBase
@reexport using Unitful
@reexport using UnPack

export Interpolations

# Common types and methods
export Layer, SubSurface, Top, Bottom
export Process, SubSurfaceProcess, BoundaryProcess, CoupledProcesses
export Coupled, Coupled2, Coupled3, Coupled4
export DiscreteEvent, ContinuousEvent, GridContinuousEvent
export VarInitializer
include("types.jl")

export Prognostic, Algebraic, Diagnostic, Var
export VarDim, OnGrid, Shape, Scalar, GridOffset, Edges, Cells
export varname, vartype, vardims, varunits, vardomain, isprognostic, isalgebraic, isflux, isdiagnostic, isongrid, dimlength
include("variables.jl")

export BCKind, Volume, FixedVolume, DiagnosticVolume, PrognosticVolume
include("traits.jl")

export initialcondition!, updatestate!, interact!, interactmaybe!, computefluxes!, resetfluxes!
export variables, processes, initializers, timestep, isactive, caninteract
export boundaryflux, boundaryvalue, criterion, criterion!, trigger!
include("methods.jl")

# Submodules

export convert_t, convert_tspan, pstrip, @pstrip, @sym_str
include("Utils/Utils.jl")
using .Utils

export DiscretizationStrategy, AutoGrid, PresetGrid, LinearSpacing, Grid, cells, edges, subgridinds, Δ, volume, area, getvar
include("Numerics/Numerics.jl")
using .Numerics

export initializer
include("initializers.jl")

include("IO/InputOutput.jl")
@reexport using .InputOutput

include("Tiles/Tiles.jl")
@reexport using .Tiles
parameters = Tiles.parameters

export ConstantBC, PeriodicBC, ConstantValue, PeriodicValue, ConstantFlux, PeriodicFlux
export volumetricfractions
include("Physics/Physics.jl")

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
