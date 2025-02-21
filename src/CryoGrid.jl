module CryoGrid

global DEBUG = haskey(ENV,"CG_DEBUG") && ENV["CG_DEBUG"] == "true"
@deprecate debug(x) debug!(x)

"""
    debug!(debug::Bool)

Enable or disable global debug mode for CryoGrid. Debug mode disables certain optimizations (e.g. loop vectorizations)
which sometimes ineterfere with the debugger. It also enables additional numerical stability checks.
"""
function debug!(debug::Bool)
    global DEBUG = debug
    # disable loop vectorization in debug mode
    Numerics.turbo(!debug)
    DEBUG && @warn "Debug mode enabled! Some performance features such as loop vectorization are now turned off by default."
    return DEBUG
end

global AUTOPARA = haskey(ENV,"CG_AUTOPARA") && ENV["CG_AUTOPARA"] == "true"

"""
    autoparam!(parameterize::Bool)

Enable or disable automatic parameterization mode. When enabled, model parameters in all layer/process types will be initialized
as free `Param`s rather than `FixedParam`s.
"""
function autoparam!(parameterize::Bool)
    global AUTOPARA = parameterize
    AUTOPARA && @warn "Automatic parameterization enabled!" || @warn "Automatic parameterization disabled!"
    return AUTOPARA
end

using Adapt
using Base: @propagate_inbounds
using ComponentArrays
using ConstructionBase
using Dates
using IfElse
using LinearAlgebra
using ModelParameters
using Reexport
using Requires

import Flatten
import Interpolations

# Re-exported third-party packages (for convenience)
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

using SciMLBase: intervals # resolve name collision

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

export BCKind
include("traits.jl")

export initialcondition!, computediagnostic!, interact!, interactmaybe!, computeprognostic!, resetfluxes!, diagnosticstep!
export variables, processes, initializers, timestep, isactive, caninteract, param
export boundaryflux, boundaryvalue, criterion, criterion!, trigger!
include("methods.jl")

# Submodules

export convert_t, convert_tspan, pstrip, @pstrip, @sym_str
include("Utils/Utils.jl")
using .Utils

export DiscretizationStrategy, AutoGrid, PresetGrid, LinearSpacing, Grid, cells, edges, subgridinds, Δ, volume, area, getvar
include("Numerics/Numerics.jl")
using .Numerics

include("default_grids.jl")

export initializer
include("initializers.jl")

include("IO/InputOutput.jl")
@reexport using .InputOutput

include("Tiles/Tiles.jl")
@reexport using .Tiles
using .Tiles: layers # resolve name collision
parameters = Tiles.parameters

export ConstantBC, PeriodicBC, ConstantValue, PeriodicValue, ConstantFlux, PeriodicFlux
export volumetricfractions
include("Physics/Physics.jl")

include("Diagnostics/Diagnostics.jl")
@reexport using .Diagnostics

# Coupling
include("coupling.jl")

# SciML problem interface
export CryoGridProblem
include("problem.jl")

# Solvers
export CGEuler, CryoGridIntegrator, CryoGridSolution
include("Solvers/Solvers.jl")

# Built-in model definitions
export SoilHeatTile, SamoylovDefault
include("models.jl")

using .InputOutput: Resource

Forcings = (
    Samoylov_ERA5_fitted_daily_1979_2020 = Resource("samoylov_era5_fitted_daily_1979-2020", ForcingFormatJSON{2}(), "https://nextcloud.awi.de/s/ScYAoHzeMzAfpjf/download/samoylov_era5_fitted_daily_1979-2020.json"),
    Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044 = Resource("Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044", ForcingFormatJSON{1}(), "https://nextcloud.awi.de/s/cbeycGQoQpXi3Ei/download/samoylov_ERA_obs_fitted_1979_2014_spinup_extended2044.json"),
    Samoylov_ERA_MkL3_CCSM4_long_term = Resource("Samoylov_ERA_MkL3_CCSM4_long_term", ForcingFormatJSON{1}(), "https://nextcloud.awi.de/s/45ax9AsTACxL25Q/download/FORCING_ULC_126_72.json"),
    Bayelva_ERA5_fitted_daily_1979_2020 = Resource("bayelva_era5_fitted_daily_1979-2020", ForcingFormatJSON{2}(), "https://nextcloud.awi.de/s/5AdbRMYKneCHgx4/download/bayelva_era5_fitted_daily_1979-2020.json")
)
Parameters = (
    # Faroux et al. doi:10.1109/IGARSS.2007.4422971
    EcoCLimMap_ULC_126_72 = Resource("EcoCLimMap_ULC_126_72", ParamsJSON{1}(), "https://nextcloud.awi.de/s/7F65JET9TzdosMD/download/PARA_ULC_126_72.json")
)

function __init__()
    # SimulationBasedInference extension
    @require SimulationBasedInference="78927d98-f421-490e-8789-96b006983a5c" begin
        using .SimulationBasedInference
        include("../ext/CryoGridSBIExt/CryoGridSBIExt.jl")
        @reexport using .CryoGridSBIExt
    end
end

end # module
