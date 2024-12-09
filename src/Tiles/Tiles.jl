module Tiles

using CryoGrid
using CryoGrid: Parameterization, DynamicParameterization, DVar
using CryoGrid.InputOutput: CryoGridParams
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

import ConstructionBase
import Flatten
import ForwardDiff
import Interpolations

export Stratigraphy, @Stratigraphy, NamedLayer
export layernames, layertypes, layers, boundaries, layername, namedlayers
export top, bottom, subsurface
include("stratigraphy.jl")

export TileState
include("state.jl")

include("tile_base.jl")

export Tile, withaxes, getstate
include("tile.jl")

export JacobianStyle, DefaultJac, TridiagJac
include("jac.jl")

precompile(CryoGridParams, (Tile,))

end