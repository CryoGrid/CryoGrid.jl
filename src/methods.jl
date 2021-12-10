"""
    variables(::Layer, ::Process)

Defines variables for a given Process on the given Layer. Implementations should return a `Tuple` of `Var`s.
"""
variables(::Layer, ::Process) = ()
"""
    callbacks(::Layer, ::Process)

Defines callbacks for a given Process on the given Layer. Implementations should return a `Tuple` or `Callback`s.
"""
callbacks(::Layer, ::Process) = ()
"""
    initialcondition!(::Layer, ::Process, state)

Defines the initial condition for a given Process on the given Layer. `initialcondition!` should write initial values into all relevant
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

Defines a boundary interaction between two processes on adjacent layers. For any interaction, the order of the arguments
follows decreasing depth, i.e. the first layer/process is always on top of the second layer/process. This ordering matters
and separate dispatches must be provided for interactions in reverse order.
"""
interact!(::Layer, ::Process, ::Layer, ::Process, state1, state2) = nothing
"""
    boundaryflux(bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, s1, s2)
    boundaryflux(s::BoundaryStyle, bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, s1, s2)

Computes the flux dH/dt at the boundary layer. Calls boundaryflux(BoundaryStyle(B),...) to allow for generic implementations by boundary condition type.
"""
boundaryflux(bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, s1, s2) = boundaryflux(BoundaryStyle(bc), bc, b, p, sub, s1, s2)
boundaryflux(s::BoundaryStyle, bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, s1, s2) = error("missing implementation of $(typeof(b)) $(typeof(s)) boundaryflux for $(typeof(bc)) + $(typeof(p)) on $(typeof(sub))")
"""
    boundaryvalue(bc::BoundaryProcess, b, p, layer, sbc, ssub)

Computes the value of the boundary condition specified by `bc` for the given layer/process combinations.
"""
boundaryvalue(bc::BoundaryProcess, b, p, layer, sbc, ssub) = error("missing implementation of boundaryvalue for $(typeof(b)) $(typeof(bc)) on $(typeof(layer)) with $(typeof(p))")
"""
    criterion(c::Callback, ::Layer, ::Process, state)

Callback criterion/condition. Should return a `Bool` for discrete callbacks and a real number for continuous callbacks.
"""
criterion(c::Callback, ::Layer, ::Process, state) = error("missing implementation of criterion for $(typeof(c))")
"""
    affect!(c::Callback, ::Layer, ::Process, state)

Callback action executed when `criterion` is met (boolean condition for discrete callbacks, zero for continuous callbacks).
"""
affect!(c::Callback, ::Layer, ::Process, state) = error("missing implementation of affect! for $(typeof(c))")
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
# default callback style impl
CallbackStyle(::C) where {C<:Callback} = CallbackStyle(C)
CallbackStyle(::Type{<:Callback}) = Discrete()
