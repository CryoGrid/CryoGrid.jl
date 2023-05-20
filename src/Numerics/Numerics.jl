module Numerics

using CryoGrid
using CryoGrid.Utils

using Base: @inbounds, @propagate_inbounds
using ComponentArrays
using DimensionalData: AbstractDimArray, DimArray, Dim, At, dims, Z
using Flatten
using IfElse
using Interpolations
using IntervalSets
using LinearAlgebra
using LoopVectorization
using Unitful
using UnPack
using StaticArrays
using StructTypes

import ConstructionBase
import ForwardDiff
import ModelParameters
import PreallocationTools as Prealloc

USE_TURBO = true
function turbo(value::Bool)
    global USE_TURBO = value
end

export Profile, ProfileKnot
include("profile.jl")

export ∇, flux!, divergence!, nonlineardiffusion!, harmonicmean!, harmonicmean
include("math.jl")

export Grid, UnitRectangle, cells, edges, subgridinds, Δ, volume, area, updategrid!
include("grid.jl")

export DiscretizationStrategy, PresetGrid, AutoGrid, makegrid, discretize
include("discretization.jl")

export DiffCache, retrieve
include("diffcache.jl")

export StateVars, getvar, getvars
include("statevars.jl")

end
