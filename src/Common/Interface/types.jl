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
"""
    Boundary = Union{Top,Bottom}

Alias that refers to the type union over both Top and Bottom layer types.
"""
const Boundary = Union{Top,Bottom}
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
    BoundaryProcess{P<:SubSurfaceProcess}

Abstract base type for boundary processes, i.e. processes that operate at the boundaries of the
subsurface. A `BoundaryProcess` represents the boundary conditions of one or more `SubSurfaceProcess`es
but may include its own diagnostic (or even prognostic) variables, if necessary. The type of `P` indicates
which `SubSurfaceProcess` is forced by the `BoundaryProcess` and may be a `Union{...}` type if the `BoundaryProcess`
provides boundary conditions to multiple `SubSurfaceProcess`es.
"""
abstract type BoundaryProcess{P<:SubSurfaceProcess} <: Process end
"""
    System{TProcs} <: Process

Represents a explicitly or implicitly coupled system of processes. `TProcs` is always a `Tuple`
of other processes.
"""
struct System{TProcs} <: Process
    processes::TProcs
    System(processes::SubSurfaceProcess...) = new{typeof(processes)}(processes)
    System(processes::BoundaryProcess...) = new{typeof(processes)}(processes)
end
"""
    Coupled{P1,P2} = System{Tuple{T1,T2}} where {T1,T2}

Represents a coupled pair of explicitly processes. Alias for `System{Tuple{P1,P2}}`.
`Coupled` provides a simple mechanism for defining new behaviors on composite processes/systems.
"""
const Coupled{P1,P2} = System{Tuple{T1,T2}} where {T1,T2}
"""
    Coupled(p1,p2)

Alias for `System(p1,p2)`.
"""
Coupled(p1::P1, p2::P2) where {P1<:Process,P2<:Process} = System(p1,p2)
# Base methods
Base.show(io::IO, ps::System{T}) where T = print(io, "$T")
@propagate_inbounds @inline Base.getindex(ps::System, i) = ps.processes[i]
# allow broadcasting of Process types
Base.Broadcast.broadcastable(p::Process) = Ref(p)

# Boundary condition trait
"""
Trait that specifies the "style" or kind of boundary condition. This can be used to write generic
implementations of `interact!` that are (relatively) agnostic to specific implementations of
`BoundaryProcess`. A good example of this can be found in `HeatConduction.boundaryflux`.
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
    AbstractParameterization

Base type for representing parameterizations.
"""
abstract type AbstractParameterization end
const Parameterization = Union{Nothing,<:AbstractParameterization}
