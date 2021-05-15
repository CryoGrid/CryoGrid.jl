# Base type for composite parameter types
abstract type Params end
# scalar broadcasting of Params types
Base.Broadcast.broadcastable(p::Params) = Ref(p)
# provide length and iteration over field names
Base.length(p::Params) = fieldcount(typeof(p))
Base.iterate(p::Params) = structiterate(p)
Base.iterate(p::Params, state) = structiterate(p,state)

export Params

# Layer base types
abstract type Layer end
abstract type SubSurface <: Layer end
struct Top <: Layer end
struct Bottom <: Layer end
const Boundary = Union{Top,Bottom}
# allow broadcasting of Layer types
Base.Broadcast.broadcastable(l::Layer) = Ref(l)

export Layer, SubSurface, Top, Bottom, Boundary

# Base types for dynamical processes
abstract type Process end
abstract type SubSurfaceProcess <: Process end
abstract type BoundaryProcess{P<:SubSurfaceProcess} <: Process end
struct Processes{TProcs} <: Process
    processes::TProcs
    Processes(processes::Process...) = new{typeof(processes)}(processes)
end
Base.show(io::IO, ps::Processes{T}) where T = print(io, "$T")
@propagate_inbounds @inline Base.getindex(ps::Processes, i) = ps.processes[i]
# allow broadcasting of Process types
Base.Broadcast.broadcastable(p::Process) = Ref(p)

export Process, SubSurfaceProcess, BoundaryProcess, Processes

# Boundary condition trait
"""
Trait that specifies the "style" or kind of boundary condition.
"""
abstract type BoundaryStyle end
struct Dirichlet <: BoundaryStyle end
struct Neumann <: BoundaryStyle end
struct Robin <: BoundaryStyle end
# Default to an error to avoid erroneously labeling a boundary condition.
BoundaryStyle(::Type{T}) where {T<:BoundaryProcess} = error("No style specified for boundary condition $T")

export BoundaryStyle, Dirichlet, Neumann, Robin
