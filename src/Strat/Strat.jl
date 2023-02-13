module Strat

import CryoGrid
import CryoGrid.InputOutput

using CryoGrid
using CryoGrid: Parameterization, DynamicParameterization
using CryoGrid.InputOutput: Forcing, CryoGridParams
using CryoGrid.Numerics
using CryoGrid.Physics
using CryoGrid.Utils

using ComponentArrays
using DataStructures: OrderedDict
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
import Interpolations
import ModelParameters: update

export TileState, LayerState
include("state.jl")

export Stratigraphy, @Stratigraphy, NamedLayer
export layernames, layertypes, layers, boundaries, layername, stratiterate
export top, bottom, subsurface
include("stratigraphy.jl")

export ConstantInitializer, InterpInitializer, initializer, init!
include("init.jl")

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
    JacobianStyle(::Type{<:Tile})

Can be overriden/extended to specify Jacobian structure for specific `Tile`s.
"""
JacobianStyle(::Type{<:Tile}) = DefaultJac()
# Auto-detect Jacobian sparsity for problems with one or more heat-only layers.
# Note: This assumes that the processes/forcings on the boundary layers do not violate the tridiagonal structure!
# Unfortunately, the Stratigraphy type signature is a bit nasty to work with :(
const HeatOnlyTile = Tile{<:Stratigraphy{N,<:Tuple{TTop,Vararg{<:Union{<:Named{<:Any,<:SubSurface{<:HeatBalance}},TBot}}}}} where {N,TTop,TBot}
JacobianStyle(::Type{<:HeatOnlyTile}) = TridiagJac()

precompile(CryoGridParams, (Tile,))

end