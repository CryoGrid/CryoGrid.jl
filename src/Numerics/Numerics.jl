module Numerics

import Base.==
import ConstructionBase
import ForwardDiff
import ModelParameters
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

struct ProfileKnot{D,T}
    depth::D
    value::T
end
Base.show(io::IO, knot::ProfileKnot) = print(io, "$(knot.depth): $(knot.value)")
struct Profile{N,TKnots}
    knots::TKnots
    Profile(::Tuple{}) = new{0,Tuple{}}(())
    Profile(knots::Tuple{Vararg{<:ProfileKnot,N}}) where N = new{N,typeof(knots)}(knots)
    Profile(pairs::Tuple{Vararg{<:Pair}}) where D = Profile(map(Base.splat(ProfileKnot), pairs))
    Profile(pairs::Pair...) = Profile(pairs)
end
Flatten.flattenable(::Type{<:ProfileKnot}, ::Type{Val{:depth}}) = false
Flatten.flattenable(::Type{ProfileKnot{D,T}}, ::Type{Val{:value}}) where {D,T<:Number} = false
Flatten.flattenable(::Type{ProfileKnot{D,T}}, ::Type{Val{:depth}}) where {D<:ModelParameters.Param,T} = true
Flatten.flattenable(::Type{ProfileKnot{D,T}}, ::Type{Val{:value}}) where {D,T<:ModelParameters.Param} = true
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

export Grid, cells, edges, subgridinds, Δ, volume, area, updategrid!
include("grid.jl")

export Var, Prognostic, Algebraic, Diagnostic, VarDim, OnGrid, Shape, Scalar
export varname, vartype, vardims, isprognostic, isalgebraic, isflux, isdiagnostic, isongrid, dimlength
include("variables.jl")

export VarStates, DiffCache, retrieve, getvar, getvars
include("varstates.jl")

include("discretize.jl")

export ConstantInitializer, InterpInitializer, initializer, init!
include("init.jl")

end
