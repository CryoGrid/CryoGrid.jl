# Layer base types
abstract type Layer end
abstract type SubSurface <: Layer end
struct Top <: Layer end
struct Bottom <: Layer end
const Boundary = Union{Top,Bottom}
# allow broadcasting of Layer types
Base.Broadcast.broadcastable(l::Layer) = Ref(l)

"""
    Process

Abstract base type for dynamical processes, e.g. heat conduction, water infiltration, etc.
"""
abstract type Process end
"""
    SubSurfaceProcess

Abstract base type for subsurface processes, i.e. processes that operate at or below the surface.
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
    Coupled{P1,P2}

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
    BoundaryStyle

Trait that specifies the "style" or kind of boundary condition.
"""
abstract type BoundaryStyle end
struct Dirichlet <: BoundaryStyle end
struct Neumann <: BoundaryStyle end
