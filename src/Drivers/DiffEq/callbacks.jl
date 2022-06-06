# callback building functions
function makecallbacks(component::StratComponent{L,P,name}, i_layer) where {L,P,name}
    cbs = []
    i_ev = 1
    for ev in CryoGrid.events(component.layer, component.processes)
        push!(cbs, _diffeqcallback(ev, Val{name}(), i_layer, i_ev))
        i_ev += 1
    end
    return Tuple(cbs)
end
function makecallbacks(tile::Tile)
    diffeq_callbacks = tuplejoin((makecallbacks(comp, i) for (i,comp) in enumerate(tile.strat))...)
    return diffeq_callbacks
end
function _criterionfunc(::Val{name}, i_layer, i_ev) where name
    function _condition(u,t,integrator)
        let tile = Tile(integrator),
            comp = tile.strat[i_layer],
            layer = comp.layer,
            process = comp.processes,
            ev = tile.events[name][i_ev],
            u = Strat.withaxes(u, tile),
            du = Strat.withaxes(get_du(integrator), tile),
            t = t;
            criterion(ev, layer, process, Strat.getstate(Val{name}(), tile, u, du, t, integrator.dt))
        end
    end
end
function _triggerfunc(::Val{name}, trig, i_layer, i_ev) where name
    _invoke_trigger!(ev, ::Nothing, layer, process, state) = trigger!(ev, layer, process, state)
    _invoke_trigger!(ev, trig::ContinuousTrigger, layer, process, state) = trigger!(ev, trig, layer, process, state)
    function _trigger!(integrator)
        let tile=Tile(integrator),
            comp = tile.strat[i_layer],
            layer = comp.layer,
            process = comp.processes,
            ev = tile.events[name][i_ev],
            u = Strat.withaxes(integrator.u, tile),
            du = Strat.withaxes(get_du(integrator), tile),
            t = integrator.t,
            state = Strat.getstate(Val{name}(), tile, u, du, t, integrator.dt);
            _invoke_trigger!(ev, trig, layer, process, state)
        end
    end
end
_diffeqcallback(::DiscreteEvent, ::Val{name}, i_layer, i_ev) where name = DiffEqCallbacks.DiscreteCallback(
    _criterionfunc(Val{name}(), i_layer, i_ev),
    _triggerfunc(Val{name}(), nothing, i_layer, i_ev);
    # todo: initialize and finalize?
)
_diffeqcallback(::ContinuousEvent, ::Val{name}, i_layer, i_ev) where name = DiffEqCallbacks.ContinuousCallback(
    _criterionfunc(Val{name}(), i_layer, i_ev),
    _triggerfunc(Val{name}(), Increasing(), i_layer, i_ev);
    affect_neg! = _triggerfunc(Val{name}(), Decreasing(), i_layer, i_ev),
    # todo: initialize and finalize?
)

"""
    (p::DiffEqCallbacks.StepsizeLimiterAffect{typeof(CryoGrid.timestep)})(integrator)

Custom implementation of `StepsizeLimiterAffect` function for `CryoGrid.timestep` that invokes the
`timestep` function with `tile,du,u,p,t` as arguments.
"""
function (p::DiffEqCallbacks.StepsizeLimiterAffect{typeof(CryoGrid.timestep)})(integrator)
    integrator.opts.dtmax = p.safety_factor*p.dtFE(integrator.sol.prob.f.f, get_du(integrator), integrator.u, integrator.p, integrator.t)
    # This part is copied from the implementation in DiffEqCallbacks
    if !integrator.opts.adaptive
        if integrator.opts.dtmax < integrator.dtcache
            integrator.dtcache = integrator.opts.dtmax
        elseif p.cached_dtcache <= integrator.opts.dtmax
            integrator.dtcache = p.cached_dtcache
        end
    end
    if p.max_step && isfinite(integrator.opts.dtmax)
        set_proposed_dt!(integrator,integrator.opts.dtmax)
        integrator.dtcache = integrator.opts.dtmax
    end
    u_modified!(integrator, false)
end