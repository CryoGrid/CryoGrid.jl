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
struct System{TProcs} <: Process
    processes::TProcs
    System(processes::Process...) = new{typeof(processes)}(processes)
end
"""
    Coupled{P1,P2}

Represents a coupled pair of processes. Alias for `System{Tuple{P1,P2}}`.
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

export Process, SubSurfaceProcess, BoundaryProcess, System, Coupled

# Boundary condition trait
"""
    BoundaryStyle

Trait that specifies the "style" or kind of boundary condition.
"""
abstract type BoundaryStyle end
struct Dirichlet <: BoundaryStyle end
struct Neumann <: BoundaryStyle end

export BoundaryStyle, Dirichlet, Neumann

"""
    JacobianStyle

Trait for indicating Jacobian sparsity of a CryoGrid ODEProblem.
"""
abstract type JacobianStyle end
struct DefaultJac <: JacobianStyle end
struct TridiagJac <: JacobianStyle end

export DefaultJac, TridiagJac

"""
    Parameterization

Base type for representing parameterizations.
"""
abstract type Parameterization end
struct Nonparametric <: Parameterization end
struct Parametric{T}
    Parametric(::T) where T = new{T}()
end

export Nonparametric, Parametric
