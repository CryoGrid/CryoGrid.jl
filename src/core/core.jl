global CRYOGRID_DEBUG = haskey(ENV,"CG_DEBUG") && ENV["CG_DEBUG"] == "true"

# Convenience constants for units
const Unit{N,D,A} = Unitful.Units{N,D,A}
const DistUnit{N} = Unitful.FreeUnits{N,Unitful.𝐋,nothing} where {N}
const DistQuantity{T,U} = Quantity{T,Unitful.𝐋,U} where {T,U<:DistUnit}
const TempUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝚯,A} where {N,A}
const TempQuantity{T,U} = Quantity{T,Unitful.𝚯,U} where {T,U<:TempUnit}
const TimeUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝐓,A} where {N,A}
const TimeQuantity{T,U} = Quantity{T,Unitful.𝐓,U} where {T,U<:TempUnit}
# Qualified name for ComponentArrays.Axis
const CAxis = ComponentArrays.Axis

export Unit, DistUnit, DistQuantity, TempUnit, TempQuantity, CAxis

# include standalone(-ish) types and methods
include("utils.jl")
include("math.jl")
include("grid.jl")
include("forcing.jl")
include("variables.jl")
include("state.jl")

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

# Model core method stubs
variables(::Layer) = ()
variables(::Layer, ::Process) = ()
initialcondition!(::Layer, state) = nothing
initialcondition!(::Layer, ::Process, state) = nothing
diagnosticstep!(l::Layer, p::Process, state) = error("no diagnostic step defined for $(typeof(l)) with $(typeof(p))")
prognosticstep!(l::Layer, p::Process, state) = error("no prognostic step defined for $(typeof(l)) with $(typeof(p))")
interact!(l1::Layer, p1::Process, l2::Layer, p2::Process, state1, state2) =
    error("no interaction defined betweeen $(typeof(l1)) with $(typeof(p1)) and $(typeof(l2)) with $(tyepof(p2))")
"""
Declares an interraction to be present between two layer/process pairs. Default rule is that all processes of the
same family (i.e. same base type name) should interact.
"""
interactrule(::Type{<:Layer}, ::Type{P1}, ::Type{<:Layer}, ::Type{P2}) where {P1<:SubSurfaceProcess,P2<:SubSurfaceProcess} =
    Base.typename(P1).wrapper == Base.typename(P2).wrapper
interactrule(::Type{L1}, ::Type{P1}, ::Type{L2}, ::Type{<:BoundaryProcess{P2}}) where {L1,L2,P1,P2} =
    interactrule(L1,P1,L2,P2)
interactrule(::Type{L1}, ::Type{<:BoundaryProcess{P1}}, ::Type{L2}, ::Type{P2}) where {L1,L2,P1,P2} =
    interactrule(L1,P1,L2,P2)

export variables, parameters, initialcondition!, diagnosticstep!, prognosticstep!, interact!, interactrule

"""
Convenience macro for setting scalar (single-element) arrays/vectors. It turns an expression of the form:
    `a.b = ...`
into
    `a.b[1] = ...`

This is primarily intended for code clarity, i.e to clearly discern scalar and non-scalar values.
"""
macro setscalar(expr)
    refexpr = expr.args[1]
    valexpr = expr.args[2]
    quote
        $(esc(refexpr))[1] = $(esc(valexpr))
    end
end

getscalar(a::AbstractArray{T,1}) where T = a[1]

export @setscalar, getscalar

# include core-dependent types/functions
include("stratigraphy.jl")
include("setup.jl")
