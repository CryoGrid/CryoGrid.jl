# callback building functions
function makecallbacks(tile::Tile)
    eventname(::Event{name}) where name = name
    isgridevent(::GridContinuousEvent) = true
    isgridevent(::Event) = false
    callbacks = []
    for (i,comp) in enumerate(tile.strat)
        events = CryoGrid.events(comp.layer, comp.process)
        layername = componentname(comp)
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
            comp = tile.strat[i_layer],
            layer = comp.layer,
            process = comp.process,
            u = Strat.withaxes(u, tile),
            du = Strat.withaxes(get_du(integrator), tile),
            t = t,
            state = Strat.getstate(Val{layername}(), tile, u, du, t, integrator.dt);
            return criterion(ev, layer, process, state)
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
    dtFE = p.safety_factor*p.dtFE(integrator.sol.prob.f.f, get_du(integrator), integrator.u, integrator.p, integrator.t)
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