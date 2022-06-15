"""
Abstract base type for all layers.
"""
abstract type Layer end
"""
    SubSurface <: Layer

Abstract base type for layers in the stratigraphy, e.g. soil, snow, pond, etc.
"""
abstract type SubSurface <: Layer end
"""
    Top <: Layer

Generic "top" layer that marks the upper boundary of the subsurface grid.
"""
struct Top <: Layer end
"""
    Bottom <: Layer

Generic "bottom" layer that marks the lower boundary of the subsurface grid.
"""
struct Bottom <: Layer end
# allow broadcasting of Layer types
Base.Broadcast.broadcastable(l::Layer) = Ref(l)

"""
Abstract base type for all dynamical processes.
"""
abstract type Process end
"""
    isless(::Process, ::Process)

Defines an ordering between processes, which can be useful for automatically enforcing
order constraints when coupling processes together.
"""
Base.isless(::Process, ::Process) = false
"""
    SubSurfaceProcess <: Process

Abstract base type for subsurface processes, i.e. processes that operate at or below the surface,
such as heat conduction, water infiltration, etc.
"""
abstract type SubSurfaceProcess <: Process end
"""
    BoundaryProcess{T<:SubSurfaceProcess}

Abstract base type for boundary processes, i.e. processes that operate at the boundaries of the
subsurface. A `BoundaryProcess` represents the boundary conditions of one or more `SubSurfaceProcess`es
but may include its own diagnostic (or even prognostic) variables, if necessary.
"""
abstract type BoundaryProcess{T<:SubSurfaceProcess} <: Process end
BoundaryProcess(::Type{T}) where {T<:SubSurfaceProcess} = BoundaryProcess{U} where {T<:U<:SubSurfaceProcess}
"""
    CoupledProcesses{TProcs} <: Process

Represents an explicitly or implicitly coupled system of processes. `TProcs` is always a `Tuple`
of other processes.
"""
struct CoupledProcesses{TProcs} <: Process
    processes::TProcs
    CoupledProcesses(processes::SubSurfaceProcess...) = CoupledProcesses(processes)
    CoupledProcesses(processes::BoundaryProcess...) = CoupledProcesses(processes)
    CoupledProcesses(processes::Tuple{Vararg{Process}}) = new{typeof(processes)}(processes)
end
"""
    Coupled2{P1,P2} = CoupledProcesses{Tuple{T1,T2}} where {T1,T2}

Represents an explicitly coupled pair of processes. Alias for `CoupledProcesses{Tuple{P1,P2}}`.
`Coupled` provides a simple mechanism for defining new behaviors on multi-processes systems.
"""
const Coupled2{P1,P2} = CoupledProcesses{Tuple{T1,T2}} where {T1,T2}
const Coupled3{P1,P2,P3} = CoupledProcesses{Tuple{T1,T2,T3}} where {T1,T2,T3}
const Coupled4{P1,P2,P3,P4} = CoupledProcesses{Tuple{T1,T2,T3,T4}} where {T1,T2,T3,T4}
"""
    Coupled(ps::Process...)

Alias for `CoupledProcesses(ps...)`.
"""
Coupled(ps::Process...) = CoupledProcesses(ps...)
# Base methods
Base.show(io::IO, ::CoupledProcesses{T}) where T = print(io, "Coupled($(join(T.parameters, " with ")))")
Base.iterate(cp::CoupledProcesses) = Base.iterate(cp.processes)
Base.iterate(cp::CoupledProcesses, state) = Base.iterate(cp.processes, state)
@propagate_inbounds Base.getindex(cp::CoupledProcesses, i) = cp.processes[i]
# allow broadcasting of Process types
Base.Broadcast.broadcastable(p::Process) = Ref(p)

# Boundary condition trait
"""
Trait that specifies the "style" or kind of boundary condition. This can be used to write generic
implementations of `interact!` that are (relatively) agnostic to specific implementations of
`BoundaryProcess`. A good example of this can be found in the `boundaryflux` method interface.
"""
abstract type BoundaryStyle end
"""
`BoundaryStyle` instance for Dirichlet boundary conditions.
"""
struct Dirichlet <: BoundaryStyle end
"""
`BoundaryStyle` instance for Neumann boundary conditions.
"""
struct Neumann <: BoundaryStyle end
"""
    Event{name}

Base type for integration "events" (i.e. callbacks). `name` should be a `Symbol`
or type which identifies the event for the purposes of dispatch on `criterion` and `trigger!`.
"""
abstract type Event{name} end
struct DiscreteEvent{name} <: Event{name}
    DiscreteEvent(name::Symbol) = new{name}()
end
struct ContinuousEvent{name} <: Event{name}
    ContinuousEvent(name::Symbol) = new{name}()
end
"""
    ContinuousTrigger

Base type for continuous trigger flags, `Increasing` and `Decreasing`, which indicate
an upcrossing of the function root (negative to positive) and a downcrossing (positive to negative)
respectively.
"""
abstract type ContinuousTrigger end
struct Increasing <: ContinuousTrigger end
struct Decreasing <: ContinuousTrigger end
"""
    Parameterization

Base type for generic parameterizations of constants or unknown quantities.
"""
abstract type Parameterization end
"""
    DynamicParameterization

Base type for dynamic parameterizations whose values may be time or state dependent.
"""
abstract type DynamicParameterization <: Parameterization end
