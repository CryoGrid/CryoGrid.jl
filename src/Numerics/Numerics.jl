module Numerics

using CryoGrid
using CryoGrid.Utils

using Adapt
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

USE_TURBO = true
function turbo(value::Bool)
    global USE_TURBO = value
end

export Profile, ProfileKnot
include("profile.jl")

export ∇, flux!, divergence!, nonlineardiffusion!, harmonicmean!, harmonicmean
include("math.jl")

export Grid, UnitRectangle, cells, edges, subgridinds, Δ, volume, area, updategrid!, currentgrid
include("grid.jl")

export DiscretizationStrategy, PresetGrid, AutoGrid, LinearSpacing, makegrid, discretize
include("discretization.jl")

export DiffCache, ArrayCache, retrieve
include("caches.jl")

export StateVars, getvar, getvars
include("statevars.jl")

end
