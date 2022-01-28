module Numerics

import Base.==
import ConstructionBase
import ForwardDiff
import PreallocationTools as Prealloc

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
using RuntimeGeneratedFunctions
using Unitful
using StructTypes
using Symbolics
using SymbolicUtils

RuntimeGeneratedFunctions.init(@__MODULE__)

export GridSpec, Edges, Cells
export Profile, ProfileKnot

abstract type AbstractDiscretization{Q,N} <: DenseArray{Q,N} end

abstract type GridSpec end
struct Edges <: GridSpec end
struct Cells <: GridSpec end

abstract type Geometry end
struct UnitVolume <: Geometry end

struct ProfileKnot{D<:DistQuantity,T}
    depth::D
    value::T
end
Base.show(io::IO, knot::ProfileKnot) = print(io, "$(knot.depth): $(knot.value)")
struct Profile{N,TKnots}
    knots::NTuple{N,TKnots}
    Profile(::Tuple{}) = new{0,Tuple{}}(())
    Profile(knots::NTuple{N,TKnots}) where {N,TKnots} = new{N,TKnots}(knots)
    Profile(pairs::NTuple{N,<:Pair{D}}) where {N,D} = Profile(map(Base.splat(ProfileKnot), pairs))
    Profile(pairs::Pair{D}...) where {D} = Profile(pairs)
end
Flatten.flattenable(::Type{<:ProfileKnot}, ::Type{Val{:depth}}) = false
Base.length(::Profile{N}) where N = N
Base.iterate(profile::Profile) = iterate(profile.knots)
Base.iterate(profile::Profile, state) = iterate(profile.knots, state)
Base.getindex(profile::Profile, itrv::Interval) = Profile(Tuple(knot for knot in profile.knots if knot.depth ∈ itrv))
Base.getindex(profile::Profile, i::Int) = profile.knots[i]
Base.getindex(profile::Profile, i) = Profile(profile.knots[i])
Base.lastindex(profile::Profile) = lastindex(profile.knots)
StructTypes.StructType(::Type{<:Profile}) = StructTypes.UnorderedStruct()

export ∇, Tabulated
include("math.jl")

export Grid, cells, edges, subgridinds, Δ, volume, area
include("grid.jl")

export Var, Prognostic, Algebraic, Diagnostic, VarDim, OnGrid, Shape, Scalar
export varname, vartype, vardims, isprognostic, isalgebraic, isflux, isdiagnostic, isongrid, dimlength
include("variables.jl")

export VarStates, DiffCache, retrieve, getvar, getvars
include("varstates.jl")

include("discretize.jl")

export ConstantInit, InterpInit, initializer, init!
include("init.jl")

end
