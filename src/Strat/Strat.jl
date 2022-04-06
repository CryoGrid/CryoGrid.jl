module Strat

import CryoGrid
import CryoGrid: variables, initialcondition!, prognosticstep!, diagnosticstep!, interact!, observe

using CryoGrid: Layer, Top, Bottom, SubSurface, Process, SubSurfaceProcess, BoundaryProcess, CoupledProcesses
using CryoGrid: Parameterization, DynamicParameterization
using CryoGrid.InputOutput: Forcing, CryoGridParams
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
using Unitful
using Reexport
using Setfield

@reexport using ModelParameters
@reexport using SimulationLogs

import Flatten
import ModelParameters: update

export Stratigraphy, @Stratigraphy
export StratComponent, componentname, copmonenttypes, components, boundaries
export top, bottom, subsurface
include("stratigraphy.jl")

export TileState, LayerState
include("state.jl")

export Tile, withaxes, getstate, parameters
include("tile.jl")

precompile(CryoGridParams, (Tile,))

end