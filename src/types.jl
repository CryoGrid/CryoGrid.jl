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
    CoupledProcesses(processes::Tuple{Vararg{Process}}) = new{typeof(processes)}(processes)
    CoupledProcesses(processes::SubSurfaceProcess...) = new{typeof(processes)}(processes)
    CoupledProcesses(processes::BoundaryProcess...) = new{typeof(processes)}(processes)
end
Base.iterate(cp::CoupledProcesses) = Base.iterate(cp.processes)
Base.iterate(cp::CoupledProcesses, state) = Base.iterate(cp.processes, state)
"""
    Coupled{P1,P2} = CoupledProcesses{Tuple{T1,T2}} where {T1,T2}

Represents an explicitly coupled pair of processes. Alias for `CoupledProcesses{Tuple{P1,P2}}`.
`Coupled` provides a simple mechanism for defining new behaviors on multi-processes systems.
"""
const Coupled{P1,P2} = CoupledProcesses{Tuple{T1,T2}} where {T1,T2}
"""
    Coupled(ps::Process...)

Alias for `CoupledProcesses(ps...)`.
"""
Coupled(ps::Process...) = CoupledProcesses(ps...)
# Base methods
Base.show(io::IO, ::CoupledProcesses{T}) where T = print(io, "Coupled($(join(T.parameters, " with ")))")
@propagate_inbounds @inline Base.getindex(ps::CoupledProcesses, i) = ps.processes[i]
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
    Callback

Base type for callback implementations.
"""
abstract type Callback end
"""
    CallbackStyle

Trait for callback types.
"""
abstract type CallbackStyle end
struct Discrete <: CallbackStyle end
struct Continuous <: CallbackStyle end
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
