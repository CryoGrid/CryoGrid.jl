module Boundaries

import CryoGrid.Interface: BoundaryStyle, variables

using CryoGrid.Interface
using CryoGrid.Numerics
using CryoGrid.Utils

using Unitful

export Constant, Periodic

struct Constant{S,T} <: BoundaryProcess
    value::T
    Constant{P,S}(value::T) where {S<:BoundaryStyle,T} = new{S,T}(value)
end

# Arguments are irrelevant for Constant, so we can just use args...
(bc::Constant)(args...) = bc.value

BoundaryStyle(::Type{<:Constant{P,S}}) where {P,S} = S()


struct Periodic{S,T} <: BoundaryProcess
    period::Float"s"
    amplitude::T
    offset::T
    Periodic{S}(period::Q, amplitude::T=one(T), offset::T=one(T)) where
        {S<:BoundaryStyle,Q<:TimeQuantity,T} =
        new{S,T}(uconvert(u"s",period) |> dustrip, amplitude, offset)
end

@inline (bc::Periodic)(t) = bc.amplitude*sin(π*(1/bc.period)*t) + bc.offset
@inline (bc::Periodic)(l1,l2,p2,s1,s2) = bc(s1.t)

BoundaryStyle(::Type{<:Periodic{S}}) where {S} = S()

include("transformed.jl")

end