module Boundaries

import CryoGrid.Interface: BoundaryStyle, variables

using CryoGrid.Interface
using CryoGrid.Numerics
using CryoGrid.Utils

using Unitful

export Constant, Periodic, Bias

"""
    struct Constant{S,T,name} <: BoundaryProcess

Constant boundary condition (of any type/unit) specified by `value`. If `name` is provided,
then `Constant` will provide a tracked parameter `name`.
"""
struct Constant{S,T,name} <: BoundaryProcess
    value::T
    Constant{S}(value::T) where {S<:BoundaryStyle,T} = new{S,T,nothing}(value)
    Constant{S}(value::T, name::Symbol) where {S<:BoundaryStyle,T} = new{S,T,name}(value)
end

# no name/parameter
(bc::Constant{S,T,nothing})(l1,p2,l2,s1,s2) where {S,T} = bc.value
# with parameter
(bc::Constant{S,T,name})(l1,p2,l2,s1,s2) where {S,T,name} = s1.params[name] |> getscalar

variables(top::Top, bc::Constant{S,T,name}) where {S,T,name} = (Parameter(name, bc.value),)
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

# Aliases
const Bias{name,T} = Constant{Dirichlet,T,name} where {name,T}
Bias(;value::T=0.0, name::Symbol=:bias) where T = Constant{Dirichlet}(value, name)

include("composite.jl")

end