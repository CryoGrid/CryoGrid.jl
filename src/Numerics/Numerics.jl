module Numerics

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
using StaticArrays
using StructTypes

import Base.==
import ConstructionBase
import ForwardDiff
import ModelParameters
import PreallocationTools as Prealloc

USE_TURBO = true
function turbo(value::Bool)
    global USE_TURBO = value
end

export Var, GridSpec, Edges, Cells, UnitVolume
include("types.jl")

export Profile, ProfileKnot
include("profile.jl")

export ∇, flux!, divergence!, nonlineardiffusion!, harmonicmean!, harmonicmean
include("math.jl")

export Grid, cells, edges, subgridinds, Δ, volume, area, updategrid!
include("grid.jl")

export Prognostic, Algebraic, Diagnostic, VarDim, OnGrid, Shape, Scalar
export varname, vartype, vardims, varunits, vardomain, isprognostic, isalgebraic, isflux, isdiagnostic, isongrid, dimlength
include("variables.jl")

export StateVars, DiffCache, retrieve, getvar, getvars
include("statevars.jl")

end
