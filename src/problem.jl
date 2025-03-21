"""
    CryoGridProblem{iip,Tu,Tt,Tp,TT,Tsv,Tsf,Tcb,Tdf,Tkw} <: SciMLBase.AbstractODEProblem{Tu,Tt,iip}

Represents a CryoGrid discretized PDE forward model configuration using the `SciMLBase`/`DiffEqBase` problem interface.
"""
struct CryoGridProblem{iip,Tu,Tt,Tp,TT,Tsv,Tsf,Tcb,Tdf,Tkw} <: SciMLBase.AbstractODEProblem{Tu,Tt,iip}
    f::TT
    u0::Tu
    tspan::NTuple{2,Tt}
    p::Tp
    callbacks::Tcb
    saveat::Tsv
    savefunc::Tsf
    isoutofdomain::Tdf
    kwargs::Tkw
    CryoGridProblem{iip}(f::TF, u0::Tu, tspan::NTuple{2,Tt}, p::Tp, cbs::Tcb, saveat::Tsv, savefunc::Tsf, iood::Tdf, kwargs::Tkw) where {iip,TF,Tu,Tt,Tp,Tsv,Tsf,Tcb,Tdf,Tkw} =
        new{iip,Tu,Tt,Tp,TF,Tsv,Tsf,Tcb,Tdf,Tkw}(f, u0, tspan, p, cbs, saveat, savefunc, iood, kwargs)
end

"""
    CryoGridProblem(tile::Tile, u0::ComponentVector, tspan::NTuple{2,DateTime}, args...;kwargs...)
"""
CryoGridProblem(tile::Tile, u0::ComponentVector, tspan::NTuple{2,DateTime}, args...;kwargs...) = CryoGridProblem(tile, u0, convert_tspan(tspan), args...;kwargs...)
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
        step_limiter=timestep,
        safety_factor=1,
        max_step=true,
        callback=nothing,
        isoutofdomain=Tiles.domain(tile),
        specialization=SciMLBase.AutoSpecialize,
        function_kwargs=(),
        prob_kwargs...
    )

Constructor for `CryoGridProblem` that automatically generates all necessary callbacks.
"""
function CryoGridProblem(
    tile::Tile,
    u0::ComponentVector,
    tspan::NTuple{2,Float64},
    p::Union{Nothing,AbstractVector}=nothing;
    diagnostic_stepsize=3600.0,
    saveat=3600.0,
    savevars=(),
    save_everystep=false,
    save_start=true,
    save_end=true,
    step_limiter=timestep,
    safety_factor=1,
    max_step=true,
    callback=nothing,
    isoutofdomain=Tiles.domain(tile),
    specialization=SciMLBase.AutoSpecialize,
    function_kwargs=(),
    prob_kwargs...
)
    getsavestate(tile::Tile, u, du) = deepcopy(Tiles.getvars(tile.state, Tiles.withaxes(u, tile), Tiles.withaxes(du, tile), savevars...))
    savefunc(u, t, integrator) = getsavestate(Tile(integrator), Tiles.withaxes(u, Tile(integrator)), get_du(integrator))
    # strip all "fixed" parameters
    tile = stripparams(FixedParam, tile)
    # retrieve variable parameters
    tilepara = length(ModelParameters.params(tile)) > 0 ? parameters(tile) : nothing
    du0 = zero(u0)
    # remove units
    tile = stripunits(tile)
    # set up saving callback
    stateproto = getsavestate(tile, u0, du0)
    savevals = SavedValues(Float64, typeof(stateproto))
    saveat = expandtstep(saveat, tspan)
    savingcallback = SavingCallback(savefunc, savevals; saveat=saveat, save_start=save_start, save_end=save_end, save_everystep=save_everystep)
    diagnostic_step_callback = PresetTimeCallback(tspan[1]:diagnostic_stepsize:tspan[end], diagnosticstep!)
    defaultcallbacks = (savingcallback, diagnostic_step_callback)
    # add step limiter to default callbacks, if defined
    if !isnothing(step_limiter)
        defaultcallbacks = (
            defaultcallbacks...,
            StepsizeLimiter(step_limiter; safety_factor, max_step)
        )
    end
    # build layer callbacks
    layercallbacks = _makecallbacks(tile)
    # add user callbacks
    usercallbacks = isnothing(callback) ? () : callback
    callbacks = CallbackSet(defaultcallbacks..., layercallbacks..., usercallbacks...)
    # note that this implicitly discards any existing saved values in the model setup's state history
    tile.data.outputs = savevals
    # build mass matrix
    mass_matrix = Numerics.build_mass_matrix(tile.state)
    # get params
    p = isnothing(p) && !isnothing(tilepara) ? ustrip.(vec(tilepara)) : p
	func = odefunction(tile, u0, p, tspan; mass_matrix, specialization, function_kwargs...)
	return CryoGridProblem{true}(func, u0, tspan, p, callbacks, saveat, getsavestate, isoutofdomain, prob_kwargs)
end

function SciMLBase.remake(
    prob::CryoGridProblem{iip};
    f=deepcopy(prob.f),
    u0=nothing,
    tspan=prob.tspan,
    p=prob.p,
    callbacks=prob.callbacks,
    saveat=prob.saveat,
    savefunc=prob.savefunc,
    isoutofdomain=prob.isoutofdomain,
    kwargs=prob.kwargs,
) where iip
    # always re-run initialcondition! with the given tspan and parameters
    _u0, du0 = initialcondition!(Tile(f), tspan, p)
    # if u0 was explicitly given, use it instead of the computed value
    if !isnothing(u0)
        # evaluate Tile on new initial state
        f(du0, u0, p, tspan[1])
    else
        u0 = _u0
    end
    return CryoGridProblem{iip}(f, u0, tspan, p, callbacks, saveat, savefunc, isoutofdomain, kwargs)
end

"""
    odefunction(setup::Tile, u0, p, tspan; kwargs...)

Constructs a SciML `ODEFunction` given the model setup, initial state u0, parameters p, and tspan.
Can (and should) be overridden by users to provide customized ODEFunction configurations for specific problem setups, e.g:
```
tile = Tile(strat,grid)
function CryoGrid.odefunction(::DefaultJac, setup::typeof(tile), u0, p, tspan)
    ...
    # make sure to return an instance of ODEFunction
end
...
prob = CryoGridProblem(tile, u0, tspan)
```

`JacobianStyle` can also be extended to create custom traits which can then be applied to compatible `Tile`s.
"""
odefunction(tile::TTile, u0, p, tspan; mass_matrix=I, specialization=SciMLBase.AutoSpecialize, kwargs...) where {TTile<:Tile} = odefunction(JacobianStyle(tile), tile, u0, p, tspan; mass_matrix, specialization, kwargs...)
odefunction(::DefaultJac, tile::Tile, u0, p, tspan; mass_matrix=I, specialization=SciMLBase.AutoSpecialize, kwargs...) = ODEFunction{true,specialization}(tile; mass_matrix, kwargs...)
function odefunction(::TridiagJac, tile::Tile, u0, p, tspan; mass_matrix=I, specialization=SciMLBase.AutoSpecialize, kwargs...)
    if :jac_prototype in keys(kwargs)
        @warn "using user specified jac_prorotype instead of tridiagonal"
        ODEFunction{true,specialization}(tile; mass_matrix, kwargs...)
    else
        N = length(u0)
        T = isnothing(p) || isa(p, SciMLBase.NullParameters) ? eltype(u0) : eltype(p)
        J = Tridiagonal(
                similar(u0, T, N-1) |> Vector,
                similar(u0, T, N) |> Vector,
                similar(u0, T, N-1) |> Vector
        )
        ODEFunction{true,specialization}(tile; jac_prototype=J, mass_matrix, kwargs...)
    end
end

# Diagnostic step callback
function diagnosticstep!(integrator::SciMLBase.DEIntegrator)
    tile = Tile(integrator)
    state = getstate(integrator)
    u_modified = diagnosticstep!(tile, state)
    # set whether or not the prognostic state u was modified
    DiffEqBase.u_modified!(integrator, u_modified)
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
    for (i,named_layer) in enumerate(namedlayers(tile.strat))
        events = CryoGrid.events(named_layer.val)
        layer_name = nameof(named_layer)
        for ev in events
            # if ev is a GridContinuousEvent, and was already added in a different layer, skip it.
            # GridContinuousEvents are defined on the whole grid/domain and so do not need to be duplicated
            # across layers.
            name = eventname(ev)
            if !isgridevent(ev) || name ∉ map(first, callbacks)
                cb = _diffeqcallback(ev, tile, Val{layer_name}(), i)
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
            u = Tiles.withaxes(u, tile),
            du = Tiles.withaxes(get_du(integrator), tile),
            t = t,
            state = getproperty(Tiles.getstate(tile, u, du, t, integrator.dt), layername);
            return criterion(ev, layer, state)
        end
    end
end

function _gridcriterionfunc(::Val{layername}, ev::Event) where layername
    function _condition(out,u,t,integrator)
        tile = Tile(integrator)
        for layer in tile.strat
            let u = Tiles.withaxes(u, tile),
                du = Tiles.withaxes(get_du(integrator), tile),
                t = t,
                state = getproperty(Tiles.getstate(tile, u, du, t, integrator.dt), layername);
                criterion!(view(out, Numerics.bounds(state.grid)), ev, layer, process, state)
            end
        end
    end
end

function _triggerfunc(::Val{layername}, ev::Event, trig::Union{Nothing,T}, i_layer::Int) where {layername,T<:ContinuousTrigger}
    _invoke_trigger!(ev, ::Nothing, layer, state) = trigger!(ev, layer, state)
    _invoke_trigger!(ev, trig::ContinuousTrigger, layer, state) = trigger!(ev, trig, layer, state)
    function _trigger!(integrator)
        let tile=Tile(integrator),
            layer = tile.strat[i_layer],
            u = Tiles.withaxes(integrator.u, tile),
            du = Tiles.withaxes(get_du(integrator), tile),
            t = integrator.t,
            state = getproperty(Tiles.getstate(tile, u, du, t, integrator.dt), layername);
            _invoke_trigger!(ev, trig, layer, state)
        end
    end
end

function _gridtriggerfunc(::Val{layername}, ev::GridContinuousEvent, grid::Grid, ::Type{T}) where {layername,T<:ContinuousTrigger}
    _invoke_trigger!(ev, ::Nothing, layer, state) = trigger!(ev, layer, state)
    _invoke_trigger!(ev, trig::ContinuousTrigger, layer, state) = trigger!(ev, trig, layer, state)
    function _trigger!(integrator, event_idx)
        tile = Tile(integrator)
        for layer in tile.strat
            u = Tiles.withaxes(integrator.u, tile)
            du = Tiles.withaxes(get_du(integrator), tile)
            t = integrator.t
            state = getproperty(Tiles.getstate(tile, u, du, t, integrator.dt), layername)
            if event_idx ∈ Numerics.bounds(state.grid)
                _invoke_trigger!(ev, T(nothing), layer, state)
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
    dtFE = p.safety_factor*p.dtFE(Tile(integrator.sol.prob.f), get_du(integrator), integrator.u, integrator.p, integrator.t)
    # unpack from ForwardDiff dual number in case autodiff is being used
    dtmax = Numerics.ForwardDiff.value(min(integrator.opts.dtmax, dtFE))
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
