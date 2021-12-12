module Strat

import CryoGrid: Layer, Top, Bottom, SubSurface, Process, SubSurfaceProcess, BoundaryProcess, CoupledProcesses
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

export TileState, LayerState
include("state.jl")

export Tile, withaxes, getstate
include("tile.jl")

export ParameterVector, LinearTrend, parameters
include("params.jl")

end