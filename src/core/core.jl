global CRYOGRID_DEBUG = haskey(ENV,"CG_DEBUG") && ENV["CG_DEBUG"] == "true"

# Convenience constants for units
const Unit{N,D,A} = Unitful.Units{N,D,A}
const DistUnit{N} = Unitful.FreeUnits{N,Unitful.𝐋,nothing} where {N}
const DistQuantity{T,U} = Quantity{T,Unitful.𝐋,U} where {T,U<:DistUnit}
# Qualified name for ComponentArrays.Axis
const CAxis = ComponentArrays.Axis

export Unit, DistUnit, DistQuantity, CAxis

# include standalone(-ish) types and methods
include("utils.jl")
include("grid.jl")
include("forcing.jl")
include("variables.jl")

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
# allow broadcasting of Layer types
Base.Broadcast.broadcastable(l::Layer) = Ref(l)

export Layer, SubSurface, Top, Bottom

# Base types for dynamical processes
abstract type Process end
abstract type SubSurfaceProcess end
abstract type BoundaryProcess{P<:SubSurfaceProcess} <: Process end
struct Processes{TProcs} <: Process
    processes::TProcs
    Processes(processes::Process...) = new{typeof(processes)}(processes)
end
# allow broadcasting of Process types
Base.Broadcast.broadcastable(p::Process) = Ref(p)

export Process, SubSurfaceProcess, BoundaryProcess, DirichletBC, NeumannBC, Processes

const State{names,T} = NamedTuple{names,T} where {names,T}
# Model core method stubs
variables(::Layer) = ()
variables(::Layer, ::Process) = ()
parameters(::Layer) = ()
parameters(::Layer, ::Process) = ()
initialcondition!(::Layer, ::State) = nothing
initialcondition!(::Layer, ::Process, ::State) = nothing
diagnosticstep!(l::Layer, p::Process, ::State) = error("no diagnostic step defined for $(typeof(l)) with $(typeof(p))")
prognosticstep!(l::Layer, p::Process, ::State) = error("no prognostic step defined for $(typeof(l)) with $(typeof(p))")
interact!(l1::Layer, p1::Process, l2::Layer, p2::Process, s1::State, s2::State) =
    error("no interaction defined betweeen $(typeof(l1)) with $(typeof(p1)) and $(typeof(l2)) with $(tyepof(p2))")
"""
Declares an interraction to be present between two layer/process pairs. Default rule is that all processes of the
same family (i.e. same base type name) should interact.
"""
interactrule(::Layer, ::P1, ::Layer, ::P2) where {P1<:SubSurfaceProcess,P2<:SubSubsurfaceProcess} =
    Base.typename(P1).wrapper == Base.typename(P2).wrapper
interactrule(::Type{P1},::Type{<:BoundaryProcess{P2}}) where {P1,P2} = interactrule(P1,P2)
interactrule(::Type{<:BoundaryProcess{P1}},::Type{P2}) where {P1,P2} = interactrule(P1,P2)

export State, variables, parameters, initialcondition!, diagnosticstep!, prognosticstep!, interact!, interactrule

"""
Convenience macro for adding dispatches. Adds CryoGrid qualifier to function names. Really just for prettiness :)
"""
macro cryogrid(func)
    @assert func.head == :function "@cryogrid may only be applied to functions"
    fname = func.args[1].args[1]
    func.args[1].args[1] = Symbol(:CryoGrid,fname)
    func
end

# include core-dependent types/functions
include("stratigraphy.jl")
include("setup.jl")
