module Sources

import CryoGrid: SubSurfaceProcess, SubSurface
import CryoGrid: diagnosticstep!, initialcondition!, interact!, prognosticstep!, variables

using ..HeatConduction
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
via `sp::S`.
"""
struct Source{P,T,S} <: SubSurfaceProcess
    term::T
    sp::S
    Source(::Type{P}, term::T, sp::S=nothing) where {P<:SubSurfaceProcess,T<:SourceTerm,S} = new{P,T,S}(term,sp)
end
ConstructionBase.constructorof(::Type{<:Source{P}}) where {P} = (term, sp) -> Source(P, term, sp)

(p::Periodic)(t) = p.amp*sin(2π*p.freq*t - p.shift) + p.level
# Heat sources
prognosticstep!(::SubSurface, s::Source{<:Heat,<:Constant}, state) = @inbounds @. state.∂H∂t += s.term.S₀
prognosticstep!(::SubSurface, src::Source{<:Heat,<:Periodic}, state) = @inbounds @. state.∂H∂t += src.term(state.t)

end