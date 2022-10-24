# Core model functions
"""
    processes(l::Layer)

Fetches the process attached to this layer, if any. Default implementation retrieves the
`proc` field on the given layer `l`. Two or more processes may be attached as a `CoupledProcess`.
"""
processes(l::Layer) = l.proc
"""
    variables(layer::Layer, process::Process)
    variables(::Layer)
    variables(::Any)

Defines variables for a given `Layer`, `Process`, or arbitrary user-defined type. Implementations should return a `Tuple` of `Var`s.
"""
variables(::Layer, process::Process) = variables(process)
variables(l::Layer) = variables(l, processes(l))
variables(::Any) = ()
"""
    basevariables(layer::Layer, process::Process)
    basevariables(::Layer)
    basevariables(::Any)

Defines "base" or common variables for a given `Layer`, `Process`, or arbitrary user-defined type. This should be used to define
mandatory or shared variables that should *not* be overridden by subtypes. As such, `basevariables` should generally only be defined
once per type hierarchy (aside from the default definitions) and only on abstract or union types.
Implementations should return a `Tuple` of `Var`s.
"""
basevariables(::Layer, process::Process) = basevariables(process)
basevariables(l::Layer) = basevariables(l, processes(l))
basevariables(::Any) = ()
"""
    initialcondition!(::Layer, state)
    initialcondition!(::Layer, ::Process, state)

Defines the initial condition for a given Process and/or Layer. `initialcondition!` should write initial values into all relevant
state variables in `state`.
"""
initialcondition!(layer::Layer, state) = initialcondition!(layer, processes(layer), state)
initialcondition!(::Layer, ::Process, state) = nothing
"""
    initialcondition!(::Layer, ::Process, ::Layer, ::Process, state1, state2)

Defines the initial condition for two processes on adjacent layers. `initialcondition!` should write initial values into all
relevant state variables in `state`.
"""
initialcondition!(layer1::Layer, layer2::Layer, state1, state2) = initialcondition!(layer1, processes(layer1), layer2, processes(layer2), state1, state2)
initialcondition!(::Layer, ::Process, ::Layer, ::Process, state1, state2) = nothing
"""
    diagnosticstep!(l::Layer, state)
    diagnosticstep!(l::Layer, p::Process, state)

Defines the diagnostic update for a Process on a given Layer.
"""
diagnosticstep!(layer::Layer, state) = diagnosticstep!(layer, processes(layer), state)
diagnosticstep!(::Layer, ::Process, state) = nothing
"""
    prognosticstep!(l::Layer, p::Process, state)

Defines the prognostic update for a Process on a given layer. Note that an instance of `prognosticstep!` must be provided
for all non-boundary (subsurface) processes/layers.
"""
prognosticstep!(layer::Layer, state) = prognosticstep!(layer, processes(layer), state)
prognosticstep!(layer::Layer, proc::Process, state) = error("no prognostic step defined for $(typeof(layer)) with $(typeof(proc))")
prognosticstep!(::Top, ::Process, state) = nothing
prognosticstep!(::Bottom, ::Process, state) = nothing
"""
    interact!(::Layer, ::Process, ::Layer, ::Process, state1, state2)

Defines a boundary interaction between two processes on adjacent layers. For any interaction, the order of the arguments
follows decreasing depth, i.e. the first layer/process is always on top of the second layer/process. This ordering matters
and separate dispatches must be provided for interactions in reverse order.
"""
interact!(layer1::Layer, layer2::Layer, state1, state2) = interact!(layer1, processes(layer1), layer2, processes(layer2), state1, state2)
interact!(::Layer, ::Process, ::Layer, ::Process, state1, state2) = nothing
"""
    timestep(::Layer, ::Process, state)

Retrieves the recommended timestep for the given `Process` defined on the given `Layer`.
The default implementation returns `Inf` which indicates no timestep restriction. The
actual chosen timestep will depend on the integrator being used and other user configuration options.
"""
timestep(layer::Layer, state) = timestep(layer, processes(layer), state)
timestep(::Layer, ::Process, state) = Inf
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
setup = Tile(stratigraphy, grid, observed=[:meanT])
# solve
...
# retrieve results
@show out.log.meanT # will be a DimArray of meanT at each timestep.
```
"""
observe(::Val{name}, layer::Layer, state) where name = observe(Val{name}(), layer, processes(layer), state)
observe(::Val{name}, ::Layer, ::Process, state) where name = nothing
# Auxiliary functions for generalized boundary implementations;
"""
    boundaryflux(bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, sbc, ssub)
    boundaryflux(s::BoundaryStyle, bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, sbc, ssub)

Computes the flux dH/dt at the boundary layer. Calls boundaryflux(BoundaryStyle(B),...) to allow for generic implementations by boundary condition type.
Note that this method uses a different argument order convention than `interact!`. This is intended to faciliate stratigraphy independent implementations
of certain boundary conditions (e.g. a simple Dirichlet boundary could be applied in the same manner to both the upper and lower boundary).
"""
boundaryflux(bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, sbc, ssub) = boundaryflux(BoundaryStyle(bc), bc, b, p, sub, sbc, ssub)
boundaryflux(s::BoundaryStyle, bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, sbc, ssub) = error("missing implementation of $(typeof(s)) $(typeof(bc)) boundaryflux on $(typeof(b)) and $(typeof(p)) on $(typeof(sub))")
"""
    boundaryvalue(bc::BoundaryProcess, lbc::Union{Top,Bottom}, proc::SubSurfaceProcess, lsub::SubSurfaceProcess, sbc, ssub)

Computes the value of the boundary condition specified by `bc` for the given layer/process combinations.
Note that this method uses a different argument order convention than `interact!`. This is intended to faciliate stratigraphy independent implementations
of certain boundary conditions (e.g. a simple Dirichlet boundary could be applied in the same manner to both the upper and lower boundary).
"""
boundaryvalue(bc::BoundaryProcess, lbc::Union{Top,Bottom}, p::SubSurfaceProcess, lsub::SubSurfaceProcess, sbc, ssub) = error("missing implementation of boundaryvalue for $(typeof(bc)) on $(typeof(lbc)) and $(typeof(p)) on $(typeof(lsub))")
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
# Events
"""
    events(::Layer, ::Process)

Defines "events" for a given Process on the given Layer. Implementations should return a `Tuple` of `Event`s.
"""
events(layer::Layer) = events(layer, processes(layer))
events(::Layer, ::Process) = ()
"""
    criterion(::Event, ::Layer, ::Process, state)

Event criterion/condition. Should return a `Bool` for discrete events. For continuous events,
this should be a real-valued function where the event is fired at the zeros/roots.
"""
criterion(ev::Union{DiscreteEvent,ContinuousEvent}, layer::Layer, state) = criterion(ev, layer, processes(layer), state)
criterion(::DiscreteEvent, ::Layer, ::Process, state) = false
criterion(::ContinuousEvent, ::Layer, ::Process, state) = Inf
"""
    criterion!(out::AbstractArray, ev::GridContinuousEvent, ::Layer, ::Process, state)

Event criterion for on-grid (i.e. multi-valued) continuous events. The condition for each grid cell should
be stored in `out`.
"""
criterion!(out::AbstractArray, ev::GridContinuousEvent, layer::Layer, state) = criterion!(out, ev, layer, processes(layer), state)
criterion!(out::AbstractArray, ::GridContinuousEvent, ::Layer, ::Process, state) = out .= Inf
"""
    trigger!(::Event, ::Layer, ::Process, state)
    trigger!(ev::ContinuousEvent, ::ContinuousTrigger, ::Layer, ::Process, state)
    trigger!(ev::GridContinuousEvent, ::ContinuousTrigger, ::Layer, ::Process, state)

Event action executed when `criterion` is met.
"""
trigger!(ev::Event, layer::Layer, state) = trigger!(ev, layer, processes(layer), state)
trigger!(::Event, ::Layer, ::Process, state) = nothing
trigger!(::ContinuousEvent, ::ContinuousTrigger, ::Layer, ::Process, state) = nothing
trigger!(::GridContinuousEvent, ::ContinuousTrigger, ::Layer, ::Process, state) = nothing
# Discretization
"""
    midpoint(::Layer, state)
    midpoint(::Layer, state, i)
    midpoint(::Layer, state, ::typeof(first))
    midpoint(::Layer, state, ::typeof(last))

Get midpoint (m) of layer or grid cell `i`.
"""
midpoint(::Layer, state) = (state.grid[end] - state.grid[1])/2
midpoint(::Layer, state, i) = cells(state.grid)[i]
midpoint(l::Layer, state, ::typeof(first)) = midpoint(l, state, 1)
midpoint(l::Layer, state, ::typeof(last)) = midpoint(l, state, lastindex(state.grid)-1)
midpoint(::Union{Top,Bottom}, state) = Inf
"""
    thickness(::Layer, state)
    thickness(::Layer, state, i)
    thickness(l::Layer, state, ::typeof(first))
    thickness(l::Layer, state, ::typeof(last))

Get thickness (m) of layer or grid cell `i`.
"""
thickness(::Layer, state) = state.grid[end] - state.grid[1]
thickness(::Layer, state, i) = Δ(state.grid)[i]
thickness(l::Layer, state, ::typeof(first)) = thickness(l, state, 1)
thickness(l::Layer, state, ::typeof(last)) = thickness(l, state, lastindex(state.grid)-1)
thickness(::Union{Top,Bottom}, state) = Inf
