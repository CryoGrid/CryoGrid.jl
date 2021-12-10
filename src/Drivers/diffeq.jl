using DiffEqBase
using DiffEqCallbacks

@reexport using OrdinaryDiffEq
@reexport using DiffEqBase: solve, init, ODEProblem, SciMLBase

"""
Specialized problem type for CryoGrid `ODEProblem`s.
"""
struct CryoGridODEProblem end

"""
    CryoGridProblem(setup::Tile, tspan::NTuple{2,Float64}, p=nothing;kwargs...)

CryoGrid specialized constructor for ODEProblem that automatically generates the initial
condition and necessary callbacks.
"""
function CryoGridProblem(
    tile::Tile,
    u0::ComponentVector,
    tspan::NTuple{2,Float64},
    p=nothing;
    saveat=3600.0,
    savevars=(),
    save_everystep=false,
    callback=nothing,
    kwargs...
)
    # workaround for bug in DiffEqCallbacks; see https://github.com/SciML/DifferentialEquations.jl/issues/326
    # we have to manually expand single-number `saveat` (i.e. time interval for saving) to a step-range.
    expandtstep(tstep::Number) = tspan[1]:tstep:tspan[end]
    expandtstep(tstep::AbstractVector) = tstep
    getsavestate(model::Tile, u, du) = deepcopy(Land.getvars(model.state, u, du, savevars...))
    savefunc(u, t, integrator) = getsavestate(integrator.f.f, Land.withaxes(u, integrator.f.f), get_du(integrator))
    pmodel = Model(tile)
    p = isnothing(p) ? dustrip.(collect(pmodel[:val])) : p
    du0 = zero(u0)
    # set up saving callback
    stateproto = getsavestate(tile, u0, du0)
    savevals = SavedValues(Float64, typeof(stateproto))
    savingcallback = SavingCallback(savefunc, savevals; saveat=expandtstep(saveat), save_everystep=save_everystep)
    layercallbacks = tuplejoin((_getcallbacks(comp) for comp in tile.strat)...)
    usercallbacks = isnothing(callback) ? () : callback
    callbacks = CallbackSet(savingcallback, layercallbacks..., usercallbacks...)
    # note that this implicitly discards any existing saved values in the model setup's state history
    tile.hist.vals = savevals
    # set up default mass matrix
    M_diag = similar(tile.state.uproto)
    M_idxmap = ComponentArrays.indexmap(getaxes(M_diag)[1])
    allvars = Flatten.flatten(tile.state.vars, Flatten.flattenable, Var)
    progvars = map(varname, filter(isprognostic, allvars))
    algvars = map(varname, filter(isalgebraic, allvars))
    for name in keys(M_idxmap)
        M_diag_var = @view M_diag[name]
        if name ∈ progvars
            M_diag_var .= one(eltype(M_diag))
        elseif name ∈ algvars
            M_diag_var .= zero(eltype(M_diag))
        end
    end
    # if no algebraic variables are present, use identity matrix
    num_algebraic = length(M_diag) - sum(M_diag)
    M = num_algebraic > 0 ? Diagonal(M_diag) : I
	func = odefunction(tile, M, u0, p, tspan; kwargs...)
	ODEProblem(func, u0, tspan, p, CryoGridODEProblem(); callback=callbacks, kwargs...)
end
"""
    CryoGridProblem(setup::Tile, tspan::NTuple{2,DateTime}, args...;kwargs...)
"""
CryoGridProblem(setup::Tile, u0::ComponentVector, tspan::NTuple{2,DateTime}, args...;kwargs...) = CryoGridProblem(setup,u0,convert_tspan(tspan),args...;kwargs...)
"""
    odefunction(setup::Tile, M, u0, p, tspan; kwargs...)

Constructs a SciML `ODEFunction` given the model setup, mass matrix M, initial state u0, parameters p, and tspan.
Can (and should) be overridden by users to provide customized ODEFunction configurations for specific problem setups, e.g:
```
tile = Tile(strat,grid)
function CryoGrid.Setup.odefunction(::DefaultJac, setup::typeof(tile), M, u0, p, tspan)
    ...
    # make sure to return an instance of ODEFunction
end
...
prob = CryoGridProblem(tile, tspan, p)
```

`JacobianStyle` can also be extended to create custom traits which can then be applied to compatible `Tile`s.
"""
odefunction(setup::TSetup, M, u0, p, tspan; kwargs...) where {TSetup<:Tile} = odefunction(JacobianStyle(TSetup), setup, M, u0, p, tspan; kwargs...)
odefunction(::DefaultJac, setup::TSetup, M, u0, p, tspan; kwargs...) where {TSetup<:Tile} = ODEFunction(setup, mass_matrix=M; kwargs...)
function odefunction(::TridiagJac, setup::Tile, M, u0, p, tspan; kwargs...)
    if :jac_prototype in keys(kwargs)
        @warn "using user specified jac_prorotype instead of tridiagonal"
        ODEFunction(setup, mass_matrix=M; kwargs...)
    else
        N = length(u0)
        J = Tridiagonal(
                similar(u0, eltype(p), N-1) |> Vector,
                similar(u0, eltype(p), N) |> Vector,
                similar(u0, eltype(p), N-1) |> Vector
        )
        ODEFunction(setup, mass_matrix=M, jac_prototype=J, kwargs...)
    end
end
"""
    getstate(layername::Symbol, integrator::SciMLBase.DEIntegrator)

Builds the state named tuple for `layername` given an initialized integrator.
"""
Land.getstate(layername::Symbol, integrator::SciMLBase.DEIntegrator) = getstate(Val{layername}(), integrator)
Land.getstate(::Val{layername}, integrator::SciMLBase.DEIntegrator) where {layername} = Land.getstate(integrator.f.f, integrator.u, get_du(integrator), integrator.t)
"""
    getvar(var::Symbol, integrator::SciMLBase.DEIntegrator)
"""
Land.getvar(var::Symbol, integrator::SciMLBase.DEIntegrator) = Land.getvar(Val{var}(), integrator.f.f, integrator.u)
"""
Constructs a `CryoGridOutput` from the given `ODESolution`.
"""
function InputOutput.CryoGridOutput(sol::TSol) where {TSol <: SciMLBase.AbstractODESolution}
    # Helper functions for mapping variables to appropriate DimArrays by grid/shape.
    withdims(::Var{name,T,<:OnGrid{Cells}}, arr, grid, ts) where {name,T} = DimArray(arr*oneunit(T), (Z(round.(typeof(1.0u"m"), cells(grid), digits=5)),Ti(ts)))
    withdims(::Var{name,T,<:OnGrid{Edges}}, arr, grid, ts) where {name,T} = DimArray(arr*oneunit(T), (Z(round.(typeof(1.0u"m"), edges(grid), digits=5)),Ti(ts)))
    withdims(::Var{name,T}, arr, zs, ts) where {name,T} = DimArray(arr*oneunit(T), (Ti(ts),))
    model = sol.prob.f.f # Tile
    ts = model.hist.vals.t # use save callback time points
    ts_datetime = Dates.epochms2datetime.(round.(ts*1000.0))
    u_all = reduce(hcat, sol.(ts))
    pax = ComponentArrays.indexmap(getaxes(model.state.uproto)[1])
    savedstates = model.hist.vals.saveval
    allvars = variables(model)
    progvars = tuplejoin(filter(isprognostic, allvars), filter(isalgebraic, allvars))
    diagvars = filter(isdiagnostic, allvars)
    fluxvars = filter(isflux, allvars)
    outputs = Dict{Symbol,Any}()
    for var in progvars
        name = varname(var)
        outputs[name] = withdims(var, u_all[pax[name],:], model.grid, ts_datetime)
    end
    for var in filter(isongrid, tuplejoin(diagvars, fluxvars))
        name = varname(var)
        states = collect(skipmissing([name ∈ keys(state) ? state[name] : missing for state in savedstates]))
        if length(states) == length(ts)
            arr = reduce(hcat, states)
            outputs[name] = withdims(var, arr, model.grid, ts_datetime)
        end
    end
    # loop over remaining (user defined) log variables
    # uservars = Dict()
    # for varname in logvarnames
    #     var_log = log[varname]
    #     if eltype(var_log) <: AbstractVector
    #         vardata = reduce(hcat, var_log)
    #         uservars[varname] = DimArray(vardata, (Z(1:size(vardata,1)),Ti(ts)))
    #     else
    #         uservars[varname] = DimArray(vardata, (Ti(ts),))
    #     end
    # end
    # if length(logvarnames) > 0
    #     nt = NamedTuple{tuple(keys(uservars)...)}(tuple(values(uservars)...))
    #     layerstates = merge(layerstates, (user=nt,))
    # end
    CryoGridOutput(ts_datetime, sol, (;outputs...))
end

"""
Evaluates the continuous solution at time `t`.
"""
(out::CryoGridOutput{<:ODESolution})(t::Real) = withaxes(out.res(t), out.res.prob.f.f)
(out::CryoGridOutput{<:ODESolution})(t::DateTime) = out(Dates.datetime2epochms(t)/1000.0)
# callback building functions
function _criterionfunc(::Val{name}, cb::Callback, layer, process) where name
    (u,t,integrator) -> let layer=layer,
        process=process,
        cb=cb,
        tile=integrator.f.f,
        u = Land.withaxes(u, tile),
        du = Land.withaxes(get_du(integrator), tile),
        t = t;
        criterion(cb, layer, process, Land.getstate(Val{name}(), tile, u, du, t))
    end
end
function _affectfunc(::Val{name}, cb::Callback, layer, process) where name
    integrator -> let layer=layer,
        process=process,
        cb=cb,
        tile=integrator.f.f,
        u = Land.withaxes(integrator.u, tile),
        du = Land.withaxes(get_du(integrator), tile),
        t = integrator.t;
        affect!(cb, layer, process, Land.getstate(Val{name}(), tile, u, du, t))
    end
end
_diffeqcallback(::Discrete, ::Val{name}, cb::Callback, layer, process) where name = DiffEqCallbacks.DiscreteCallback(
    _criterionfunc(Val{name}(), cb, layer, process),
    _affectfunc(Val{name}(), cb, layer, process),
    # todo: initialize and finalize?
)
_diffeqcallback(::Continuous, ::Val{name}, cb::Callback, layer, process) where name = DiffEqCallbacks.ContinuousCallback(
    _criterionfunc(Val{name}(), cb, layer, process),
    _affectfunc(Val{name}(), cb, layer, process),
    # todo: initialize and finalize?
)
_getcallbacks(component::StratComponent{L,P,name}) where {L,P,name} = Tuple(_diffeqcallback(CallbackStyle(callback), Val{name}(), callback, component.layer, proc) for proc in component.process for callback in callbacks(component.layer, proc))
