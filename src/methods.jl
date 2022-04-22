# Core model functions
"""
    variables(::Layer)
    variables(::Process)
    variables(layer::Layer, process::Process)

Defines variables for a given Process and/or Layer. Implementations should return a `Tuple` of `Var`s.
"""
variables(::Layer) = ()
variables(::Process) = ()
variables(layer::Layer, process::Process) = tuple(variables(layer)..., variables(process)...)
"""
    initialcondition!(::Layer, state)
    initialcondition!(::Process, state) = nothing
    initialcondition!(::Layer, ::Process, state)

Defines the initial condition for a given Process and/or Layer. `initialcondition!` should write initial values into all relevant
state variables in `state`.
"""
initialcondition!(::Layer, state) = nothing
initialcondition!(::Process, state) = nothing
initialcondition!(layer::Layer, process::Process, state) = nothing
"""
    initialcondition!(::Layer, ::Process, ::Layer, ::Process, state1, state2)

Defines the initial condition for two processes on adjacent layers. `initialcondition!` should write initial values into all
relevant state variables in `state`.
"""
initialcondition!(::Layer, ::Process, ::Layer, ::Process, state1, state2) = nothing
"""
    diagnosticstep!(l::Layer, p::Process, state)
    diagnosticstep!(l::Layer, state)

Defines the diagnostic update for a Process on a given Layer.
"""
diagnosticstep!(l::Layer, p::Process, state) = nothing
diagnosticstep!(l::Layer, state) = nothing
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
# Auxiliary functions for generalized boundary implementations;
# Note that these methods use a different argument order convention than `interact!`. This is intended to
# faciliate stratigraphy independent implementations of certain boundary conditions (e.g. a simple Dirichlet
# boundary could be applied in the same manner to both the upper and lower boundary).
"""
    boundaryflux(bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, sbc, ssub)
    boundaryflux(s::BoundaryStyle, bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, sbc, ssub)

Computes the flux dH/dt at the boundary layer. Calls boundaryflux(BoundaryStyle(B),...) to allow for generic implementations by boundary condition type.
"""
boundaryflux(bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, sbc, ssub) = boundaryflux(BoundaryStyle(bc), bc, b, p, sub, sbc, ssub)
boundaryflux(s::BoundaryStyle, bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, sbc, ssub) = error("missing implementation of $(typeof(s)) $(typeof(bc)) boundaryflux on $(typeof(b)) and $(typeof(p)) on $(typeof(sub))")
"""
    boundaryvalue(bc::BoundaryProcess, lbc::Union{Top,Bottom}, proc::SubSurfaceProcess, lsub::SubSurfaceProcess, sbc, ssub)

Computes the value of the boundary condition specified by `bc` for the given layer/process combinations.
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
# Callbacks
"""
    callbacks(::Layer, ::Process)

Defines callbacks for a given Process on the given Layer. Implementations should return a `Tuple` or `Callback`s.
"""
callbacks(::Layer, ::Process) = ()
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
    CallbackStyle(::C)
    CallbackStyle(::Type{<:Callback})

Trait implementation that defines the "style" or type of the given callback as any subtype of `CallbackStyle`,
for example `Discrete` or `Continuous`.
"""
CallbackStyle(::C) where {C<:Callback} = CallbackStyle(C)
CallbackStyle(::Type{<:Callback}) = Discrete()
# Discretization
"""
    midpoints(::Layer, state)

Get midpoint(s) of layer or all grid cells within the layer.
"""
midpoints(::Layer, state) = cells(state.grid)
"""
    thickness(::Layer, state)

Get thickness (m) of layer or all grid cells within the layer.
"""
thickness(::Layer, state) = Δ(state.grid)
# Composition
"""
    volumetricfractions(::Layer, ::Process, state)
    volumetricfractions(::Layer, ::Process, state, i)

Get the volumetric fractions of each constituent in the volume (at grid cell `i`, if specificed).
"""
volumetricfractions(::Layer, ::Process, state) = ()
volumetricfractions(::Layer, ::Process, state, i) = ()
