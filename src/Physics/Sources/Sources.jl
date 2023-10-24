module Sources

import CryoGrid: SubSurfaceProcess, SubSurface
import CryoGrid: computediagnostic!, initialcondition!, interact!, computefluxes!, variables

using ..Heat
using CryoGrid.Numerics
using CryoGrid.Utils

using ConstructionBase
using ModelParameters
using Unitful

export Source, SourceTerm

"""
Abstract base type for source terms which define a `Source` process.
"""
abstract type SourceTerm end
"""
    Constant <: SourceTerm

Parametric source term that is constant through time and space.
"""
Base.@kwdef struct Constant{S} <: SourceTerm
    S₀::S = 0.0
end
"""
    Periodic <: SourceTerm

Parametric source term that is periodic through time and constant through space.
"""
Base.@kwdef struct Periodic{A,F,S,L} <: SourceTerm
    amp::A = 1.0
    freq::F = 1.0/(3600*24)
    shift::S = 0.0
    level::L = 0.0
end
"""
    Source{P,T,S} <: SubSurfaceProcess

Generic "source" process to provide additive fluxes for one or more prognostic variables in process `P`.
The behavior is governed by the `SourceTerm`, `T`, and further user specialization can be implemented
via `aux::S`.
"""
struct Source{P,T,S} <: SubSurfaceProcess
    term::T
    aux::S
    Source(::Type{P}, term::T, aux::S=nothing) where {P<:SubSurfaceProcess,T<:SourceTerm,S} = new{P,T,S}(term,aux)
end
ConstructionBase.constructorof(::Type{<:Source{P}}) where {P} = (term, aux) -> Source(P, term, aux)

(p::Periodic)(t) = p.amp*sin(2π*p.freq*t - p.shift) + p.level
# HeatBalance sources
computefluxes!(::SubSurface, s::Source{<:HeatBalance,<:Constant}, state) = @inbounds @. state.dH += s.term.S₀
computefluxes!(::SubSurface, src::Source{<:HeatBalance,<:Periodic}, state) = @inbounds @. state.dH += src.term(state.t)

end