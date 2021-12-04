module Land

import CryoGrid: Layer, Top, Bottom, SubSurface, Process, SubSurfaceProcess, BoundaryProcess, CompoundProcess
import CryoGrid: variables, initialcondition!, prognosticstep!, diagnosticstep!, interact!, observe

using CryoGrid.Layers
using CryoGrid.Numerics
using CryoGrid.Physics
using CryoGrid.Utils

using ComponentArrays
using ConstructionBase
using DataStructures: OrderedDict
using Dates
using DimensionalData
using IfElse
using IntervalSets
using Lazy: @>>, @>, groupby
using LinearAlgebra
using Parameters
using Unitful
using Reexport
using Setfield

import Flatten

@reexport using ModelParameters
@reexport using SimulationLogs

export Stratigraphy, @Stratigraphy
export StratComponent, componentname, copmonenttypes, components, boundaries
export top, bottom, subsurface
include("stratigraphy.jl")

export LandModelState, LayerState
include("state.jl")

export LandModel, withaxes, getstate
include("model.jl")

export ParameterVector, LinearTrend, parameters
include("params.jl")

end