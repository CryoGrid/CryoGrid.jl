"""
Abstract base type for all dynamical processes.
"""
abstract type Process end

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

"""
    CoupledProcesses{TProcs} <: Process

Represents an explicitly or implicitly coupled system of processes. `TProcs` is always a `Tuple`
of other processes.
"""
struct CoupledProcesses{TProcs} <: Process
    processes::TProcs
    CoupledProcesses(processes::Process...) = new{typeof(processes)}(processes)
end

"""
    Coupled2{P1,P2} = CoupledProcesses{Tuple{T1,T2}} where {T1,T2}

Type alias for coupled processes, i.e. `CoupledProcesses{Tuple{P1,P2}}`.
`Coupled` provides a simple mechanism for defining new behaviors on multi-processes systems.
"""
const Coupled2{P1,P2} = CoupledProcesses{Tuple{T1,T2}} where {T1,T2}
const Coupled3{P1,P2,P3} = CoupledProcesses{Tuple{T1,T2,T3}} where {T1,T2,T3}
const Coupled4{P1,P2,P3,P4} = CoupledProcesses{Tuple{T1,T2,T3,T4}} where {T1,T2,T3,T4}

Coupled() = CoupledProcesses()

"""
    Coupled(ps::Process...)

Constructs a composite/coupled process from one or more processes. Alias for `CoupledProcesses(ps...)`.
"""
Coupled(ps::Process...) = CoupledProcesses(ps...)

"""
    Coupled(types::Type{<:Process}...)

Convenince method which constructs a `CoupledProcesses` type corresponding to each type in `types`, e.g:

```
Coupled(SnowMassBalance, HeatBalance) = CoupledProcesses{Tuple{T1,T2}} where {T1<:SnowMassBalance, T2<:HeatBalance}
```

also equivalent to `Coupled2{<:SnowMassBalance,<:HeatBalance}`.
"""
@generated function Coupled(types::Type{<:Process}...)
    typenames = map(i -> Symbol(:T,i), 1:length(types))
    :(CoupledProcesses{Tuple{$(typenames...)}} where {$(map(i -> :($(typenames[i]) <: $(types[i].parameters[1])), 1:length(types))...)})
end

# Base methods
Base.show(io::IO, ::CoupledProcesses{T}) where T = print(io, "Coupled($(join(T.parameters, " with ")))")
Base.iterate(cp::CoupledProcesses) = Base.iterate(cp.processes)
Base.iterate(cp::CoupledProcesses, state) = Base.iterate(cp.processes, state)
Base.length(cp::CoupledProcesses) = length(cp.processes)
Base.firstindex(cp::CoupledProcesses) = firstindex(cp.processes)
Base.lastindex(cp::CoupledProcesses) = lastindex(cp.processes)
@propagate_inbounds Base.getindex(cp::CoupledProcesses, i) = cp.processes[i]
Base.Broadcast.broadcastable(p::Process) = Ref(p) # allow broadcasting of Process types

# Layers
"""
    Layer

Abstract base type for all layers.
"""
abstract type Layer end

"""
    SubSurface <: Layer

Abstract base type for layers in the stratigraphy, e.g. soil, snow, pond, etc.
"""
abstract type SubSurface <: Layer end

"""
    Top{TProc} <: Layer

Generic "top" layer that marks the upper boundary of the subsurface grid.
"""
struct Top{TProc} <: Layer
    proc::TProc
    Top(proc::CoupledProcesses{<:Tuple{Vararg{BoundaryProcess}}}) = new{typeof(proc)}(proc)
    Top(proc::BoundaryProcess) = new{typeof(proc)}(proc)
    # convenience constructor that automatically couples the processes
    Top(procs::BoundaryProcess...) = Top(Coupled(procs...))
end

"""
    Bottom{TProc} <: Layer

Generic "bottom" layer that marks the lower boundary of the subsurface grid.
"""
struct Bottom{TProc} <: Layer
    proc::TProc
    Bottom(proc::CoupledProcesses{<:Tuple{Vararg{BoundaryProcess}}}) = new{typeof(proc)}(proc)
    Bottom(proc::BoundaryProcess) = new{typeof(proc)}(proc)
    # convenience constructor that automatically couples the processes
    Bottom(procs::BoundaryProcess...) = Bottom(Coupled(procs...))
end
# allow broadcasting of Layer types
Base.Broadcast.broadcastable(l::Layer) = Ref(l)

# Events
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

struct GridContinuousEvent{name} <: Event{name}
    GridContinuousEvent(name::Symbol) = new{name}()
end

"""
    ContinuousTrigger

Base type for continuous trigger flags, `Increasing` and `Decreasing`, which indicate
an upcrossing of the function root (negative to positive) and a downcrossing (positive to negative)
respectively. Both subtypes have a field `idx` which, for `GridContinuousEvent` is set to the grid
cell index for which the event was triggered.
"""
abstract type ContinuousTrigger end

"""
Trigger for the criterion function crossing zero from negative to positive.
"""

struct Increasing <: ContinuousTrigger
    idx::Union{Nothing,Int}
end

"""
Trigger for the criterion function crossing zero from positive to negative.
"""
struct Decreasing <: ContinuousTrigger
    idx::Union{Nothing,Int}
end

"""
    Parameterization

Base type for generic parameterizations of processes or components.
"""
abstract type Parameterization end

"""
    DynamicParameterization

Base type for dynamic parameterizations whose values may be time or state dependent.
"""
abstract type DynamicParameterization <: Parameterization end
