module Boundaries

import CryoGrid: BoundaryProcess, BoundaryStyle, Dirichlet, Neumann, Top
import CryoGrid: variables

using CryoGrid.Numerics
using CryoGrid.Utils

using ConstructionBase
using Dates
using Flatten
using Interpolations
using ModelParameters
using Parameters
using TimeSeries
using Unitful

import Flatten: flattenable

export Constant, Periodic, Bias, BoundaryEffect, Damping
export Forcing, TimeSeriesForcing, ForcingData

"""
    struct Constant{S,T} <: BoundaryProcess

Constant boundary condition (of any type/unit) specified by `value`.
"""
struct Constant{S,T} <: BoundaryProcess
    value::T
    Constant(::Type{S}, value::T) where {S<:BoundaryStyle,T} = new{S,T}(value)
end
ConstructionBase.constructorof(::Type{<:Constant{S}}) where {S} = value -> Constant(S,value)
(bc::Constant{S,T})(l1,p2,l2,s1,s2) where {S,T} = bc.value

BoundaryStyle(::Type{<:Constant{S}}) where {S} = S()

"""
    struct Periodic{S,T} <: BoundaryProcess

Periodic boundary condition (of any type/unit) specified by `period`, `amplitude`, and `phaseshift`.
"""
struct Periodic{S,T} <: BoundaryProcess
    period::Float"s"
    amplitude::T
    phaseshift::T
    Periodic{S}(period::Q, amplitude::T=one(T), phaseshift::T=one(T)) where
        {S<:BoundaryStyle,Q<:TimeQuantity,T} =
        new{S,T}(uconvert(u"s",period) |> dustrip, amplitude, phaseshift)
end

@inline (bc::Periodic)(t) = bc.amplitude*sin(π*(1/bc.period)*t + bc.phaseshift)
@inline (bc::Periodic)(l1,l2,p2,s1,s2) = bc(s1.t)

BoundaryStyle(::Type{<:Periodic{S}}) where {S} = S()

@with_kw struct Bias{P} <: BoundaryProcess
    bias::P = Param(0.0)
end
(bc::Bias)(l1,p2,l2,s1,s2) = bc.bias

BoundaryStyle(::Type{<:Bias}) = Dirichlet()

abstract type BoundaryEffect end
"""
    Damping{D,K} <: BoundaryEffect

Generic implementation of bulk conductive damping at the boundary.
"""
@with_kw struct Damping{D,K} <: BoundaryEffect
    depth::D = t -> 0.0 # function of t -> damping depth; defaults to zero function
    k::K = Param(1.0, bounds=(0.0,Inf)) # conductivity of medium
end
function (damp::Damping)(u_top, u_sub, k_sub, Δsub, t)
    let d_med = damp.depth(t), # length of medium
        d_sub = Δsub / 2, # half length of upper grid cell
        d = d_med + d_sub, # total flux distance (from grid center to "surface")
        k_med = damp.k, # conductivity of damping medium (assumed uniform)
        k_sub = k_sub, # conductivity at grid center
        k = Numerics.harmonicmean(k_med, k_sub, d_med, d_sub);
        -k*(u_sub - u_top)/d/Δsub # flux per unit volume
    end
end

include("composed.jl")
include("forcing.jl")

end