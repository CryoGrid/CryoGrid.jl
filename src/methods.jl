# Core model functions
"""
    processes(l::Layer)

Fetches the process(es) attached to this layer, if any. Returned value must be of type `Process`. If the layer has more than one process,
they should be combined together with `Coupled(procs...)`.
"""
processes(l::Layer) = error("not implemented for layer of type $(typeof(l)); maybe you forgot to add a method CryoGrid.processes(::$(nameof(typeof(l)))) = ...?")
processes(top::Top) = top.proc
processes(bot::Bottom) = bot.proc

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
    initialcondition!(::Layer, state)
    initialcondition!(::Layer, ::Process, state)
    initialcondition!(::VarInitializer, ::Layer, ::Process, state)

Defines the initial condition for a given `Layer` and possibly an `initializer`.
`initialcondition!` should compute initial values into all relevant state variables in `state`.
"""
initialcondition!(layer::Layer, state) = initialcondition!(layer, processes(layer), state)
initialcondition!(initializer::VarInitializer, layer::Layer, state) = initialcondition!(initializer, layer, processes(layer), state)
initialcondition!(::Layer, ::Process, state) = nothing
initialcondition!(::VarInitializer, ::Layer, ::Process, state) = nothing

"""
    initialcondition!(layer1::Layer, layer2::Layer, state1, state2)
    initialcondition!(::Layer, ::Process, ::Layer, ::Process, state1, state2)

Defines the initial condition for two processes on adjacent layers. `initialcondition!` should write initial values into all
relevant state variables in `state`.
"""
initialcondition!(layer1::Layer, layer2::Layer, state1, state2) = initialcondition!(layer1, processes(layer1), layer2, processes(layer2), state1, state2)
initialcondition!(initializer::VarInitializer, layer1::Layer, layer2::Layer, state1, state2) = initialcondition!(initializer, layer1, processes(layer1), layer2, processes(layer2), state1, state2)
initialcondition!(::Layer, ::Process, ::Layer, ::Process, state1, state2) = nothing

"""
    updatestate!(l::Layer, state)
    updatestate!(l::Layer, p::Process, state)

Updates all diagnostic/non-flux state variables for the given `Layer` based on the current prognostic state.
"""
updatestate!(layer::Layer, state) = updatestate!(layer, processes(layer), state)
updatestate!(::Layer, ::Process, state) = nothing

"""
    interact!(::Layer, ::Process, ::Layer, ::Process, state1, state2)

Defines a boundary interaction between two processes on adjacent layers. For any interaction, the order of the arguments
follows decreasing depth, i.e. the first layer/process is always on top of the second layer/process. This ordering matters
and separate dispatches must be provided for interactions in reverse order.
"""
interact!(layer1::Layer, layer2::Layer, state1, state2) = interact!(layer1, processes(layer1), layer2, processes(layer2), state1, state2)
interact!(::Layer, ::Process, ::Layer, ::Process, state1, state2) = nothing

"""
    computefluxes!(l::Layer, p::Process, state)

Calculates all internal fluxes for a given layer. Note that an instance of `computefluxes!` must be provided
for all non-boundary (subsurface) processes/layers.
"""
computefluxes!(layer::Layer, state) = computefluxes!(layer, processes(layer), state)
computefluxes!(layer::Layer, proc::Process, state) = error("no prognostic step defined for $(typeof(layer)) with $(typeof(proc))")
computefluxes!(::Top, ::BoundaryProcess, state) = nothing
computefluxes!(::Bottom, ::BoundaryProcess, state) = nothing

"""
    caninteract(layer1::Layer, layer2::Layer, state1, state2)
    caninteract(l1::Layer, ::Process, l2::Layer, ::Process, state1, state2)

Returns `true` if and only if the given layer/process types are able to interact based on the current state.
Defaults to checking whether both layers are currently active. This behavior should be overridden by subtypes where necessary.
"""
caninteract(layer1::Layer, layer2::Layer, state1, state2) = caninteract(layer1, processes(layer1), layer2, processes(layer2), state1, state2)
caninteract(l1::Layer, ::Process, l2::Layer, ::Process, state1, state2) = isactive(l1, state1) && isactive(l2, state2)

"""
    interactmaybe!(layer1::Layer, layer2::Layer, state1, state2)
    interactmaybe!(layer1::Layer, p1::Process, layer2::Layer, p2::Process, state1, state2)

Conditionally invokes `interact!` if and only if `caninteract` is true.
"""
interactmaybe!(layer1::Layer, layer2::Layer, state1, state2) = interactmaybe!(layer1, processes(layer1), layer2, processes(layer2), state1, state2)
function interactmaybe!(layer1::Layer, p1::Process, layer2::Layer, p2::Process, state1, state2)
    if caninteract(layer1, p1, layer2, p2, state1, state2)
        interact!(layer1, p1, layer2, p2, state1, state2)
        return true
    end
    return false
end

"""
    resetfluxes!(layer::Layer, state)
    resetfluxes!(layer::Layer, ::Process, state)

Resets all flux variables for the given layer/process to zero.
"""
resetfluxes!(layer::Layer, state) = resetfluxes!(layer, processes(layer), state)
resetfluxes!(layer::Layer, proc::Process, state) = nothing

"""
    isactive(::Layer, state)

Returns a boolean whether or not this layer is currently active in the stratigraphy and should interact with other layers.
Note that `updatestate!` and `computefluxes!` are always invoked regardless of the current state of `isactive`.
The default implementation of `isactive` always returns `true`.
"""
isactive(::Layer, state) = true

"""
    timestep(::Layer, ::Process, state)

Retrieves the recommended timestep for the given `Process` defined on the given `Layer`.
The default implementation returns `Inf` which indicates no timestep restriction. The
actual chosen timestep will depend on the integrator being used and other user configuration options.
"""
timestep(layer::Layer, state) = timestep(layer, processes(layer), state)
timestep(::Layer, ::Process, state) = Inf

"""
    parameterize(x::T) where {T}
    parameterize(x::Unitful.AbstractQuantity; props...)
    parameterize(p::Param; ignored...)

Recursively wraps `x` or nested numeric quantities in `x` with `Param` to mark them as parameters.
If `x` is already a `Param` type, `x` will be returned as-is.
If `x` is a numeric type, `x` will be wrapped in `Param` with associated properties `props`.
If `x` is a struct type, `x` will be recursively unpacked and `parameterize` called on each field.
"""
parameterize(x::Number; props...) = Param(x; props...)
parameterize(x::Unitful.AbstractQuantity; props...) = Param(ustrip(x); untis=unit(x), props...)
parameterize(p::Param; ignored...) = p
parameterize(f::Function; ignored...) = f
function parameterize(x::T; props...) where {T}
    # get field names of T, if available
    T_fieldnames = isabstracttype(T) ? () : fieldnames(T)
    # invoke parameterize on all fields
    new_fields = map(T_fieldnames) do fieldname
        fieldvalue = getfield(x, fieldname)
        parameterize(fieldvalue)
    end
    ctor = ConstructionBase.constructorof(T)
    return ctor(new_fields...)
end

"""
    initializers(::Layer)
    initializers(::Layer, ::Process)

Optional method that can be used to provide default initializers for state variables that will be run before user provided ones.
"""
initializers(layer::Layer) = initializers(layer, processes(layer))
initializers(::Layer, ::Process) = ()

# Auxiliary functions for generalized boundary implementations;
"""
    boundaryflux(bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, sbc, ssub)
    boundaryflux(s::BCKind, bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, sbc, ssub)

Computes the flux dH/dt at the boundary layer. Calls boundaryflux(BCKind(B),...) to allow for generic implementations by boundary condition type.
Note that this method uses a different argument order convention than `interact!`. This is intended to faciliate stratigraphy independent implementations
of certain boundary conditions (e.g. a simple Dirichlet boundary could be applied in the same manner to both the upper and lower boundary).
"""
boundaryflux(bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, sbc, ssub) = boundaryflux(BCKind(bc), bc, b, p, sub, sbc, ssub)
boundaryflux(::Neumann, bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, sbc, ssub) = boundaryvalue(bc, sbc)
boundaryflux(s::BCKind, bc::BoundaryProcess, b::Union{Top,Bottom}, p::SubSurfaceProcess, sub::SubSurface, sbc, ssub) = error("missing implementation of $(typeof(s)) $(typeof(bc)) boundaryflux on $(typeof(b)) and $(typeof(p)) on $(typeof(sub))")

"""
    boundaryvalue(bc::BoundaryProcess, state)

Computes the value of the boundary condition specified by `bc` for the given layer/process combinations.
Note that this method uses a different argument order convention than `interact!`. This is intended to faciliate stratigraphy independent implementations
of certain boundary conditions (e.g. a simple Dirichlet boundary could be applied in the same manner to both the upper and lower boundary).
"""
boundaryvalue(bc::BoundaryProcess, state) = error("missing implementation of boundaryvalue for $(typeof(bc))")

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
