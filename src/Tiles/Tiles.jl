module Tiles

import CryoGrid
import CryoGrid.InputOutput

using CryoGrid
using CryoGrid: Parameterization, DynamicParameterization
using CryoGrid.InputOutput: Forcing, CryoGridParams
using CryoGrid.Numerics
using CryoGrid.Utils

using Adapt
using ComponentArrays
using Dates
using DimensionalData
using IfElse
using IntervalSets
using LinearAlgebra
using Unitful
using Reexport
using Setfield

@reexport using ModelParameters
@reexport using SimulationLogs

import ConstructionBase
import Flatten
import ForwardDiff
import Interpolations

export TileState, LayerState
include("state.jl")

export Stratigraphy, @Stratigraphy, NamedLayer
export layernames, layertypes, layers, boundaries, layername, namedlayers
export top, bottom, subsurface
include("stratigraphy.jl")

export Tile, withaxes, getstate
include("tile.jl")

export JacobianStyle, DefaultJac, TridiagJac, HeatOnlyTile

"""
    JacobianStyle

Trait for indicating Jacobian sparsity of a CryoGrid ODEProblem.
"""
abstract type JacobianStyle end
struct DefaultJac <: JacobianStyle end
struct TridiagJac <: JacobianStyle end
"""
    JacobianStyle(::Tile)

Can be overriden/extended to specify Jacobian structure for specific `Tile`s.
"""
function JacobianStyle(tile::Tile)
    prognostic_grid_vars = filter(isongrid, filter(isprognostic, CryoGrid.variables(tile)))
    if length(prognostic_grid_vars) == 1
        # Auto-detect Jacobian sparsity for problems with one on-grid prognostic variable
        return TridiagJac()
    else
        return DefaultJac()
    end
end

precompile(CryoGridParams, (Tile,))

end