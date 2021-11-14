module Numerics

import Base.==
import ExprTools

using CryoGrid.Utils

using Base: @inbounds, @propagate_inbounds
using ConstructionBase
using DimensionalData: DimArray, Dim, At, dims, Z
using Flatten
using IfElse
using Interpolations: Interpolations, Gridded, Linear, Flat, Line, interpolate, extrapolate
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

abstract type AbstractDiscretization{Q,N} <: DenseArray{Q,N} end

abstract type GridSpec end
struct Edges <: GridSpec end
struct Cells <: GridSpec end

abstract type Geometry end
struct UnitVolume <: Geometry end

export ∇
include("math.jl")

export Grid, cells, edges, indexmap, subgrid, Δ, volume, area
include("grid.jl")

export Profile, profile2array, interpolateprofile!
include("profile.jl")

export Var, Prognostic, Algebraic, Diagnostic, VarDim, OnGrid, Shape, Scalar
export varname, vartype, isprognostic, isalgebraic, isdiagnostic
include("variables.jl")

"""
    discretize([::Type{A}], ::T,  ::Var) where {T,N,D<:AbstractDiscretization{T,N},A<:AbstractArray{T,N}}

Produces a discretization of the given variable based on `T` and array type `A`.
"""
discretize(::Type{A}, ::D, ::Var) where {Q,T,N,D<:AbstractDiscretization{Q,N},A<:AbstractArray{T,N}} = error("missing discretize implementation for $D")
discretize(d::AbstractDiscretization{Q,N}, var::Var) where {Q,N} = discretize(Array{vartype(var),N}, d, var)
# grid discretizations
discretize(::Type{A}, grid::Grid, var::Var) where {A<:AbstractVector} = similar(A{vartype(var)}, dimlength(var.dim, grid))
function discretize(::Type{A}, grid::Grid, pvars::Union{<:Prognostic,<:Algebraic}...) where {A<:AbstractVector}
    # separate into grid and non-grid vars
    gridvars = filter(v -> isa(vardims(v), OnGrid), pvars)
    pointvars = filter(v -> !isa(vardims(v), OnGrid), pvars)
    # get lengths
    gridvar_ns = map(dimlength ∘ vardims, gridvars)
    pointvar_ns = map(dimlength ∘ vardims, pointvars)
    Ng, Np = sum(gridvar_ns), sum(pointvar_ns)
    # build axis indices;
    # non-grid prognostic variables get collected at the top of the vector, in the order provided
    pointvar_ax = (;(varname(p) => i:(i+n) for (p,n,i) in zip(pointvars, pointvar_ns, cumsum(vcat([1],pointvar_ns[1:end-1]))))...)
    # grid variables get interlaced throughout the rest of the vector; i.e. for variable i, its grid points are:
    # i:k:kn where k is the number of grid variables and n is the length of the grid.
    gridvar_ax = (;(varname(p) => st:length(gridvars):Ng for (p,st) in zip(gridvars, Np:(Np+length(gridvars))))...)
    # allocate component array; assumes all variables have (and should!) have the same type
    u = similar(A{vartype(first(pvars))}, Ng+Np)
    ComponentVector(u, (Axis(merge(pointvar_ax, gridvar_ax)),))
end

end
