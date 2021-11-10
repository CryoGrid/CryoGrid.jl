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

Defines the initial condition for two processes on adjacent layers. `initialcondition!` should write initial values into all
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
"""
    observe(::Val{name}, ::Layer, ::Process, state1)

Called at the end of each step. Can be used by the user to add additional observables via `@log` without affecting the
model implementation. As such, this function should **not** be used in implementations, but only by users and code which
monitors/consumes CryoGrid model outputs.

Example:
```julia
observe(::Val{:meanT}, ::SubSurface, ::Heat, state) = @log meanT = mean(state.T)
# build model
...
setup = LandModel(stratigraphy, grid, observed=[:meanT])
# solve
...
# retrieve results
@show out.log.meanT # will be a DimArray of meanT at each timestep.
```
"""
observe(::Val{name}, ::Layer, ::Process, state) where name = nothing
"""
    BoundaryStyle(::Type{T})

Can be overriden by `BoundaryProcess` types to indicate the type of boundary condition, e.g:

```
BoundaryStyle(::Type{BP}) = Dirichlet()
```

where `BP` is a `BoundaryProcess` that provides the boundary conditions.
"""
BoundaryStyle(::Type{BP}) where {BP<:BoundaryProcess} = error("No style specified for boundary process $BP")
BoundaryStyle(bc::BoundaryProcess) = BoundaryStyle(typeof(bc))
