module Sources

import CryoGrid: SubSurfaceProcess, SubSurface
import CryoGrid: diagnosticstep!, initialcondition!, interact!, prognosticstep!, variables

using ..HeatConduction
using CryoGrid.Numerics
using CryoGrid.Utils

export Source, SourceTerm

"""
Abstract base type for source terms which define a `Source` process.
"""
abstract type SourceTerm end
"""
    Constant{name} <: SourceTerm

Parametric source term with `name` that is constant through time and space.
"""
struct Constant{name} <: SourceTerm
    Constant(name::Symbol=:S₀) = new{name}()
end
"""
    Periodic{a,f,s} <: SourceTerm

Parametric source term with `name` that is periodic through time and constant through space.
`a`, `f`, and `s` are the parameter names for amplitude, frequency, and phase shift respecitvely,
each of which will be prefixed by `name`.
"""
struct Periodic{T,a,f,s} <: SourceTerm
    center::T
    Periodic(name::Symbol=:S, center::T=0.0) where T = new{T,Symbol(name,:_amp),Symbol(name,:_freq),Symbol(name,:_shift)}(center)
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

variables(::SubSurface, ::Source{P,Constant{name}}) where {P,name} = (
    Parameter(name, 0.0),
)
variables(::SubSurface, ::Source{P,Periodic{T,a,f,s}}) where {P,T,a,f,s} = (
    Parameter(a, 1.0),
    Parameter(f, 1.0),
    Parameter(s, 0.0),
)
(p::Periodic)(a,f,t₀,c,t) = a*sin(2π*f*t - t₀) + c
# Heat sources
function prognosticstep!(::SubSurface, ::Source{<:Heat,Constant{name}}, state) where {name}
    @inbounds @. state.dH += state.params[name]
end
function prognosticstep!(::SubSurface, src::Source{<:Heat,Periodic{T,a,f,s}}, state) where {T,a,f,s}
    let p = src.term,
        t = state.t,
        c = p.center
        A = state.params[a]
        ω = state.params[f]
        t₀ = state.params[s];
        @inbounds @. state.dH += p(A,ω,t₀,c,t)
    end
end

end