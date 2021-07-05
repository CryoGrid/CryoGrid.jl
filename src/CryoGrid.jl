module CryoGrid

global CRYOGRID_DEBUG = haskey(ENV,"CG_DEBUG") && ENV["CG_DEBUG"] == "true"

using Reexport

# Independent submodules
include("common/Common.jl")
include("io/InputOutput.jl")

@reexport using .Common

# Top-level method stubs
"""
    variables(::Layer)

Defines variables for a given Layer type. Implementations should return a `Tuple` of `CryoGrid.Common.Var`s.
"""
variables(::Layer) = ()
"""
    variables(::Layer, ::Process)

Defines variables for a given Process on Layer. Implementations should return a `Tuple` of `CryoGrid.Common.Var`s.
"""
variables(::Layer, ::Process) = ()
"""
    initialcondition!(::Layer, state)

Defines the initial condition for a given Layer. `initialcondition!` should write initial values into all relevant
state variables in `state`.
"""
initialcondition!(::Layer, state) = nothing
"""
    initialcondition!(::Layer, ::Process, state)

Defines the initial condition for a given Process on Layer. `initialcondition!` should write initial values into all relevant
state variables in `state`.
"""
initialcondition!(::Layer, ::Process, state) = nothing
"""
    initialcondition!(::Layer, ::Process, ::Layer, ::Process, state1, state2)

Defines the initial condition for two Processes on adjacent layers. `initialcondition!` should write initial values into all
relevant state variables in `state`.
"""
initialcondition!(::Layer, ::Process, ::Layer, ::Process, state1, state2) = nothing
"""
    diagnosticstep!(l::Layer, p::Process, state)

Defines the diagnostic update for a Process on a given Layer.
"""
diagnosticstep!(l::Layer, p::Process, state) = nothing
"""
    prognosticstep!(l::Layer, p::Process, state)

Defines the prognostic update for a Process on a given layer. Note that an instance of `prognosticstep!` must be provided
for all non-boundary (subsurface) processes/layers.
"""
prognosticstep!(l::Layer, p::Process, state) = error("no prognostic step defined for $(typeof(l)) with $(typeof(p))")
"""
    interact!(::Layer, ::Process, ::Layer, ::Process, state1, state2)

Defines a boundary interaction between two processes on adjacent layers.
"""
interact!(::Layer, ::Process, ::Layer, ::Process, state1, state2) = nothing

# We need to include the Setup module here because the JacobianStyle method definition depends on CryoGridSetup
include("setup/Setup.jl")

@reexport using .Setup

"""
    BoundaryStyle(::Type{T})

Can be overriden by `BoundaryProcess` types to indicate the type of boundary condition, e.g:

```
BoundaryStyle(::Type{BP}) = Dirichlet()
```

where `BP` is a `BoundaryProcess` that provides the boundary conditions.
"""
BoundaryStyle(::Type{T}) where {T<:BoundaryProcess} = error("No style specified for boundary process $T")
"""
    JacobianStyle(::Type{<:CryoGridSetup})

Can be overriden/extended to specify Jacobian structure for specific `CryoGridSetup`s.
"""
JacobianStyle(::Type{<:CryoGridSetup}) = DefaultJac()

export variables, initialcondition!, diagnosticstep!, prognosticstep!, interact!, BoundaryStyle, JacobianStyle

"""
    parameters(setup::CryoGridSetup; unconstrained=false)

Helper function to obtain the parameters of a `CryoGridSetup`. If `unconstrained=true`, the parameters will be mapped
to an unconstrained space first via `unconstrain`. Otherwise, they will be left as their default/initialized values.
"""
parameters(setup::CryoGridSetup; unconstrained=false) = unconstrained ? unconstrain(copy(setup.pproto), setup) : copy(setup.pproto)

# Dependent submodules
include("layers/Layers.jl")
include("processes/Processes.jl")
include("callbacks/Callbacks.jl")

# Re-export submodules
@reexport using .InputOutput: loadforcings, InputSpec, JsonSpec
@reexport using .Layers
@reexport using .Processes
@reexport using .Callbacks

# Include Models submodule last to allow dependence on other submodules.
include("models/Models.jl")

const CryoGridModels = Models
export CryoGridModels

end # module
