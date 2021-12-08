module Boundaries

import CryoGrid: BoundaryProcess, BoundaryStyle, Dirichlet, Neumann, Top
import CryoGrid: variables, boundaryvalue

using CryoGrid.Numerics
using CryoGrid.Utils

using Base: @propagate_inbounds
using ConstructionBase
using Dates
using Flatten
using Interpolations
using ModelParameters
using Parameters
using TimeSeries
using Unitful

import Flatten: flattenable

export Constant, Periodic, Bias
export BoundaryEffect, Damping
include("effects.jl")
export Forcing, TimeSeriesForcing, ForcingData
include("forcing.jl")

"""
    struct Constant{S,T} <: BoundaryProcess

Constant boundary condition (of any type/unit) specified by `value`.
"""
struct Constant{S,T} <: BoundaryProcess
    value::T
    Constant(::Type{S}, value::T) where {S<:BoundaryStyle,T} = new{S,T}(value)
end
ConstructionBase.constructorof(::Type{<:Constant{S}}) where {S} = value -> Constant(S,value)
boundaryvalue(bc::Constant{S,T},l1,p2,l2,s1,s2) where {S,T} = bc.value

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

@inline boundaryvalue(bc::Periodic,l1,p2,l2,s1,s2) = bc.amplitude*sin(π*(1/bc.period)*t + bc.phaseshift)

BoundaryStyle(::Type{<:Periodic{S}}) where {S} = S()

@with_kw struct Bias{P} <: BoundaryProcess
    bias::P = Param(0.0)
end
@inline boundaryvalue(bc::Bias,l1,p2,l2,s1,s2) = bc.bias

BoundaryStyle(::Type{<:Bias}) = Dirichlet()

"""
    struct CombinedBoundaryProcess{B1,B2,F,S} <: BoundaryProcess

Represents a composition of two boundary processes, `B1` and `B2`, via an operator `F`.
A typical use case is combining `Constant` with a forcing-driven boundary process to
scale or shift the forcing.
"""
struct CombinedBoundaryProcess{B1,B2,F,S} <: BoundaryProcess
    op::F
    bc1::B1
    bc2::B2
    function CombinedBoundaryProcess(op::F, bc1::B1, bc2::B2) where {F,B1<:BoundaryProcess,B2<:BoundaryProcess}
        @assert BoundaryStyle(bc1) == BoundaryStyle(bc2) "boundary condition styles (e.g. Dirichlet vs Neumann) must match"
        new{B1,B2,F,typeof(BoundaryStyle(bc1))}(op,bc1,bc2)
    end
end
@inline boundaryvalue(cbc::CombinedBoundaryProcess,l1,p2,l2,s1,s2) = cbc.op(cbc.bc1(l1,l2,p2,s1,s2), cbc.bc2(l1,l2,p2,s1,s2))
variables(top::Top, cbc::CombinedBoundaryProcess) = tuplejoin(variables(top, cbc.bc1), variables(top, cbc.bc2))
BoundaryStyle(::Type{CombinedBoundaryProcess{B1,B2,F,S}}) where {F,B1,B2,S} = S()
# Overload arithmetic operators on boundary processes.
Base.:+(bc1::BoundaryProcess, bc2::BoundaryProcess) = CombinedBoundaryProcess(+, bc1, bc2)
Base.:-(bc1::BoundaryProcess, bc2::BoundaryProcess) = CombinedBoundaryProcess(-, bc1, bc2)
Base.:*(bc1::BoundaryProcess, bc2::BoundaryProcess) = CombinedBoundaryProcess(*, bc1, bc2)
Base.:/(bc1::BoundaryProcess, bc2::BoundaryProcess) = CombinedBoundaryProcess(/, bc1, bc2)

end
