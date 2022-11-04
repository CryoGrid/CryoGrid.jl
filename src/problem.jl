using DiffEqCallbacks, SciMLBase
import DiffEqBase, DiffEqCallbacks

"""
    CryoGridProblem{Tu,Tt,Tp,TT,Tcb,Tdf,Tkw} <: SciMLBase.AbstractODEProblem{Tu,Tt,true}

Represents a CryoGrid discretized PDE forward model configuration using the `SciMLBase`/`DiffEqBase` problem interface.
"""
struct CryoGridProblem{Tu,Tt,Tp,TT,Tsf,Tcb,Tdf,Tkw} <: SciMLBase.AbstractODEProblem{Tu,Tt,true}
    f::TT
    u0::Tu
    tspan::NTuple{2,Tt}
    p::Tp
    callbacks::Tcb
    savefunc::Tsf
    isoutofdomain::Tdf
    kwargs::Tkw
end
"""
    CryoGridProblem(
        tile::Tile,
        u0::ComponentVector,
        tspan::NTuple{2,Float64},
        p=nothing;
        saveat=3600.0,
        savevars=(),
        save_everystep=false,
        save_start=true,
        save_end=true,
        step_limiter=CryoGrid.timestep,
        safety_factor=1,
        max_step=true,
        callback=nothing,
        isoutofdomain=Strat.domain(tile),
        specialization=SciMLBase.AutoSpecialize,
        function_kwargs=(),
        prob_kwargs...
    )

Constructor for `CryoGridProblem` that automatically generates necessary callbacks for saving
diagnostic state variables
"""
function CryoGridProblem(
    tile::Tile,
    u0::ComponentVector,
    tspan::NTuple{2,Float64},
    p=nothing;
    saveat=3600.0,
    savevars=(),
    save_everystep=false,
    save_start=true,
    save_end=true,
    step_limiter=CryoGrid.timestep,
    safety_factor=1,
    max_step=true,
    callback=nothing,
    isoutofdomain=Strat.domain(tile),
    specialization=SciMLBase.AutoSpecialize,
    function_kwargs=(),
    prob_kwargs...
)
    # workaround for bug in DiffEqCallbacks; see https://github.com/SciML/DifferentialEquations.jl/issues/326
    # we have to manually expand single-number `saveat` (i.e. time interval for saving) to a step-range.
    expandtstep(tstep::Number) = tspan[1]:tstep:tspan[end]
    expandtstep(tstep::AbstractVector) = tstep
    getsavestate(tile::Tile, u, du) = deepcopy(Strat.getvars(tile.state, Strat.withaxes(u, tile), Strat.withaxes(du, tile), savevars...))
    savefunc(u, t, integrator) = getsavestate(Tile(integrator), Strat.withaxes(u, Tile(integrator)), get_du(integrator))
    tile, p = if isnothing(p) && isempty(ModelParameters.params(tile))
        tile, nothing
    else
        model_tile = Model(tile)
        p = isnothing(p) ? collect(model_tile[:val]) : p
        model_tile[:val] = p
        parent(model_tile), p
    end
    du0 = zero(u0)
    # remove units
    tile = stripunits(tile)
    # set up saving callback
    stateproto = getsavestate(tile, u0, du0)
    savevals = SavedValues(Float64, typeof(stateproto))
    savingcallback = SavingCallback(savefunc, savevals; saveat=expandtstep(saveat), save_start=save_start, save_end=save_end, save_everystep=save_everystep)
    # add step limiter to default callbacks, if defined
    defaultcallbacks = isnothing(step_limiter) ? (savingcallback,) : (savingcallback, StepsizeLimiter(step_limiter; safety_factor, max_step))
    # build layer callbacks
    layercallbacks = _makecallbacks(tile)
    # add user callbacks
    usercallbacks = isnothing(callback) ? () : callback
    callbacks = CallbackSet(defaultcallbacks..., layercallbacks..., usercallbacks...)
    # note that this implicitly discards any existing saved values in the model setup's state history
    tile.hist.vals = savevals
    M = Numerics.build_mass_matrix(tile.state)
	func = odefunction(tile, u0, p, tspan; mass_matrix=M, specialization, function_kwargs...)
	return CryoGridProblem(func, u0, tspan, p, callbacks, getsavestate, isoutofdomain, prob_kwargs)
end
"""
    CryoGridProblem(tile::Tile, u0::ComponentVector, tspan::NTuple{2,DateTime}, args...;kwargs...)
"""
CryoGridProblem(tile::Tile, u0::ComponentVector, tspan::NTuple{2,DateTime}, args...;kwargs...) = CryoGridProblem(tile, u0, convert_tspan(tspan), args...;kwargs...)
"""
    odefunction(setup::Tile, u0, p, tspan; kwargs...)

Constructs a SciML `ODEFunction` given the model setup, initial state u0, parameters p, and tspan.
Can (and should) be overridden by users to provide customized ODEFunction configurations for specific problem setups, e.g:
```
tile = Tile(strat,grid)
function CryoGrid.Setup.odefunction(::DefaultJac, setup::typeof(tile), u0, p, tspan)
    ...
    # make sure to return an instance of ODEFunction
end
...
prob = CryoGridProblem(tile, tspan, p)
```

`JacobianStyle` can also be extended to create custom traits which can then be applied to compatible `Tile`s.
"""
odefunction(tile::TTile, u0, p, tspan; mass_matrix=I, specialization=SciMLBase.AutoSpecialize, kwargs...) where {TTile<:Tile} = odefunction(JacobianStyle(TTile), tile, u0, p, tspan; mass_matrix, specialization, kwargs...)
odefunction(::DefaultJac, tile::Tile, u0, p, tspan; mass_matrix=I, specialization=SciMLBase.AutoSpecialize, kwargs...) = ODEFunction{true,specialization}(tile; mass_matrix, kwargs...)
function odefunction(::TridiagJac, tile::Tile, u0, p, tspan; mass_matrix=I, specialization=SciMLBase.AutoSpecialize, kwargs...)
    if :jac_prototype in keys(kwargs)
        @warn "using user specified jac_prorotype instead of tridiagonal"
        ODEFunction{true,specialization}(tile; mass_matrix, kwargs...)
    else
        N = length(u0)
        J = Tridiagonal(
                similar(u0, eltype(p), N-1) |> Vector,
                similar(u0, eltype(p), N) |> Vector,
                similar(u0, eltype(p), N-1) |> Vector
        )
        ODEFunction{true,specialization}(tile; jac_prototype=J, mass_matrix, kwargs...)
    end
end
# overrides to make SciML problem interface work
SciMLBase.ODEProblem(prob::CryoGridProblem) = ODEProblem(prob.f, prob.u0, prob.tspan, prob.p; callback=prob.callbacks, isoutofdomain=prob.isoutofdomain, prob.kwargs...)
DiffEqBase.get_concrete_problem(prob::CryoGridProblem, isadapt; kwargs...) = prob
# callback building functions
function _makecallbacks(tile::Tile)
    eventname(::Event{name}) where name = name
    isgridevent(::GridContinuousEvent) = true
    isgridevent(::Event) = false
    callbacks = []
    for (i,named_layer) in enumerate(tile.strat)
        events = CryoGrid.events(named_layer.obj)
        layername = Strat.layername(named_layer)
        for ev in events
            # if ev is a GridContinuousEvent, and was already added in a different layer, skip it.
            # GridContinuousEvents are defined on the whole grid/domain and so do not need to be duplicated
            # across layers.
            name = eventname(ev)
            if !isgridevent(ev) || name ∉ map(first, callbacks)
                cb = _diffeqcallback(ev, tile, Val{layername}(), i)
                push!(callbacks, eventname(ev) => cb)
            end
        end
    end
    return map(last, callbacks)
end
function _criterionfunc(::Val{layername}, ev::Event, i_layer::Int) where layername
    function _condition(u,t,integrator)
        let tile = Tile(integrator),
            layer = tile.strat[i_layer],
            u = Strat.withaxes(u, tile),
            du = Strat.withaxes(get_du(integrator), tile),
            t = t,
            state = Strat.getstate(Val{layername}(), tile, u, du, t, integrator.dt);
            return criterion(ev, layer.obj, state)
        end
    end
end
function _gridcriterionfunc(::Val{layername}, ev::Event) where layername
    function _condition(out,u,t,integrator)
        tile = Tile(integrator)
        for comp in tile.strat
            let layer = comp.layer,
                process = comp.process,
                u = Strat.withaxes(u, tile),
                du = Strat.withaxes(get_du(integrator), tile),
                t = t,
                state = Strat.getstate(Val{layername}(), tile, u, du, t, integrator.dt);
                criterion!(view(out, Numerics.bounds(state.grid)), ev, layer, process, state)
            end
        end
    end
end
function _triggerfunc(::Val{layername}, ev::Event, trig::Union{Nothing,T}, i_layer::Int) where {layername,T<:ContinuousTrigger}
    _invoke_trigger!(ev, ::Nothing, layer, process, state) = trigger!(ev, layer, process, state)
    _invoke_trigger!(ev, trig::ContinuousTrigger, layer, process, state) = trigger!(ev, trig, layer, process, state)
    function _trigger!(integrator)
        let tile=Tile(integrator),
            comp = tile.strat[i_layer],
            layer = comp.layer,
            process = comp.process,
            u = Strat.withaxes(integrator.u, tile),
            du = Strat.withaxes(get_du(integrator), tile),
            t = integrator.t,
            state = Strat.getstate(Val{layername}(), tile, u, du, t, integrator.dt);
            _invoke_trigger!(ev, trig, layer, process, state)
        end
    end
end
function _gridtriggerfunc(::Val{layername}, ev::GridContinuousEvent, grid::Grid, ::Type{T}) where {layername,T<:ContinuousTrigger}
    _invoke_trigger!(ev, ::Nothing, layer, process, state) = trigger!(ev, layer, process, state)
    _invoke_trigger!(ev, trig::ContinuousTrigger, layer, process, state) = trigger!(ev, trig, layer, process, state)
    function _trigger!(integrator, event_idx)
        tile = Tile(integrator)
        for comp in tile.strat
            u = Strat.withaxes(integrator.u, tile)
            du = Strat.withaxes(get_du(integrator), tile)
            t = integrator.t
            state = Strat.getstate(Val{layername}(), tile, u, du, t, integrator.dt)
            if event_idx ∈ Numerics.bounds(state.grid)
                _invoke_trigger!(ev, T(nothing), comp.layer, comp.process, state)
                break
            end
        end
    end
end
_diffeqcallback(ev::DiscreteEvent, ::Tile, ::Val{layername}, i_layer::Int) where {layername} = DiffEqCallbacks.DiscreteCallback(
    _criterionfunc(Val{layername}(), ev, i_layer),
    _triggerfunc(Val{layername}(), ev, nothing, i_layer),
    # todo: initialize and finalize?
)
_diffeqcallback(ev::ContinuousEvent, ::Tile, ::Val{layername}, i_layer::Int) where layername = DiffEqCallbacks.ContinuousCallback(
    _criterionfunc(Val{layername}(), ev, i_layer),
    _triggerfunc(Val{layername}(), ev, Increasing(nothing), i_layer),
    _triggerfunc(Val{layername}(), ev, Decreasing(nothing), i_layer),
    # todo: initialize and finalize?
)
_diffeqcallback(ev::GridContinuousEvent, tile::Tile, ::Val{layername}, ::Int) where layername = DiffEqCallbacks.VectorContinuousCallback(
    _gridcriterionfunc(Val{layername}(), ev),
    _gridtriggerfunc(Val{layername}(), ev, tile.grid, Increasing),
    _gridtriggerfunc(Val{layername}(), ev, tile.grid, Decreasing),
    length(cells(tile.grid)),
    # todo: initialize and finalize?
)

"""
    (p::DiffEqCallbacks.StepsizeLimiterAffect{typeof(CryoGrid.timestep)})(integrator)

Custom implementation of `StepsizeLimiterAffect` function for `CryoGrid.timestep` that invokes the
`timestep` function with `tile,du,u,p,t` as arguments.
"""
function (p::DiffEqCallbacks.StepsizeLimiterAffect{typeof(CryoGrid.timestep)})(integrator)
    dtFE = p.safety_factor*p.dtFE(Strat.Tile(integrator.sol.prob.f), get_du(integrator), integrator.u, integrator.p, integrator.t)
    dtmax = min(integrator.opts.dtmax, dtFE)
    # This part is copied from the implementation in DiffEqCallbacks
    if !integrator.opts.adaptive
        if dtmax < integrator.dtcache
            integrator.dtcache = dtmax
        elseif p.cached_dtcache <= dtmax
            integrator.dtcache = p.cached_dtcache
        end
    end
    if p.max_step && isfinite(dtmax)
        set_proposed_dt!(integrator, dtmax)
        integrator.dtcache = dtmax
    end
    u_modified!(integrator, false)
end