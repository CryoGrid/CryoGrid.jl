"""
Driver module SciML diffeq solvers.
"""
module DiffEq

using CryoGrid.Drivers
using CryoGrid: Strat, SubSurface, CoupledProcesses, Callback, CallbackStyle, Discrete, Continuous
using CryoGrid.InputOutput
using CryoGrid.Numerics
using CryoGrid.Physics: Heat
using CryoGrid.Utils

import CryoGrid: variables, callbacks, criterion, affect!
import CryoGrid.Strat: Tile, Stratigraphy, StratComponent

using ComponentArrays
using Dates
using DimensionalData
using Flatten
using IfElse
using IntervalSets
using ModelParameters
using LinearAlgebra
using LinearSolve
using Reexport
using Unitful

using DiffEqBase
using DiffEqBase.SciMLBase
using DiffEqCallbacks

@reexport using OrdinaryDiffEq
@reexport using DiffEqBase: solve, init, ODEProblem, SciMLBase

export CryoGridProblem

export CFLStepLimiter
include("steplimiters.jl")

export TDMASolver
include("solvers.jl")

"""
Specialized problem type for CryoGrid `ODEProblem`s.
"""
struct CryoGridODEProblem end

function Strat.Tile(integrator::SciMLBase.DEIntegrator)
    tile = integrator.sol.prob.f.f
    return Strat.updateparams(tile, Strat.withaxes(integrator.u, tile), integrator.p, integrator.t)
end

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
    save_start=true,
    save_end=true,
    callback=nothing,
    kwargs...
)
    # workaround for bug in DiffEqCallbacks; see https://github.com/SciML/DifferentialEquations.jl/issues/326
    # we have to manually expand single-number `saveat` (i.e. time interval for saving) to a step-range.
    expandtstep(tstep::Number) = tspan[1]:tstep:tspan[end]
    expandtstep(tstep::AbstractVector) = tstep
    getsavestate(tile::Tile, u, du) = deepcopy(Strat.getvars(tile.state, Strat.withaxes(u, tile), Strat.withaxes(du, tile), savevars...))
    savefunc(u, t, integrator) = getsavestate(Tile(integrator), Strat.withaxes(u, Tile(integrator)), get_du(integrator))
    # remove units
    tile = stripunits(tile)
    # collect parameters
    p = isnothing(p) ? dustrip.(collect(Model(tile)[:val])) : collect(p)
    du0 = zero(u0)
    # set up saving callback
    stateproto = getsavestate(tile, u0, du0)
    savevals = SavedValues(Float64, typeof(stateproto))
    savingcallback = SavingCallback(savefunc, savevals; saveat=expandtstep(saveat), save_start=save_start, save_end=save_end, save_everystep=save_everystep)
    layercallbacks = _makecallbacks(tile)
    usercallbacks = isnothing(callback) ? () : callback
    callbacks = CallbackSet(savingcallback, layercallbacks..., usercallbacks...)
    # note that this implicitly discards any existing saved values in the model setup's state history
    tile.hist.vals = savevals
    # set up default mass matrix, M:
    # M⋅∂u∂t = f(u)
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
    getstate(integrator::SciMLBase.DEIntegrator)
    getstate(layername::Symbol, integrator::SciMLBase.DEIntegrator)

Builds the state named tuple for `layername` given an initialized integrator.
"""
Strat.getstate(integrator::SciMLBase.DEIntegrator) = Strat.getstate(Tile(integrator), integrator.u, get_du(integrator), integrator.t)
Strat.getstate(layername::Symbol, integrator::SciMLBase.DEIntegrator) = Strat.getstate(Val{layername}(), integrator)
Strat.getstate(::Val{layername}, integrator::SciMLBase.DEIntegrator) where {layername} = Strat.getstate(Val{layername}(), Tile(integrator), integrator.u, get_du(integrator), integrator.t)
"""
    getvar(var::Symbol, integrator::SciMLBase.DEIntegrator)
"""
Numerics.getvar(var::Symbol, integrator::SciMLBase.DEIntegrator) = Numerics.getvar(Val{var}(), Tile(integrator), integrator.u)
"""
    CryoGridOutput(sol::TSol, tspan::NTuple{2,Float64}=(-Inf,Inf)) where {TSol<:SciMLBase.AbstractODESolution}

Constructs a `CryoGridOutput` from the given `ODESolution`. Optional argument `tspan` restricts the time span of the output.
"""
InputOutput.CryoGridOutput(sol::TSol, tspan::NTuple{2,DateTime}) where {TSol<:SciMLBase.AbstractODESolution} = CryoGridOutput(sol, convert_tspan(tspan))
function InputOutput.CryoGridOutput(sol::TSol, tspan::NTuple{2,Float64}=(-Inf,Inf)) where {TSol<:SciMLBase.AbstractODESolution}
    # Helper functions for mapping variables to appropriate DimArrays by grid/shape.
    withdims(var::Var{name,<:OnGrid{Cells}}, arr, grid, ts) where {name} = DimArray(arr*one(vartype(var))*varunits(var), (Z(round.(typeof(1.0u"m"), cells(grid), digits=5)),Ti(ts)))
    withdims(var::Var{name,<:OnGrid{Edges}}, arr, grid, ts) where {name} = DimArray(arr*one(vartype(var))*varunits(var), (Z(round.(typeof(1.0u"m"), edges(grid), digits=5)),Ti(ts)))
    withdims(var::Var{name}, arr, zs, ts) where {name} = DimArray(arr*one(vartype(var))*varunits(var), (Ti(ts),))
    save_interval = ClosedInterval(tspan...)
    model = sol.prob.f.f # Tile
    ts = model.hist.vals.t # use save callback time points
    t_mask = map(∈(save_interval), ts) # indices within t interval
    u_all = reduce(hcat, sol.(ts[t_mask])) # build prognostic state from continuous solution
    pax = ComponentArrays.indexmap(getaxes(model.state.uproto)[1])
    # get saved diagnostic states and timestamps only in given interval
    savedstates = model.hist.vals.saveval[t_mask]
    ts_datetime = Dates.epochms2datetime.(round.(ts[t_mask]*1000.0))
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
        if length(states) == length(ts_datetime)
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
function _criterionfunc(::Val{name}, i_layer, i_proc, i_cb) where name
    function _condition(u,t,integrator)
        let tile = Tile(integrator),
            comp = tile.strat[i_layer],
            layer = comp.layer,
            process = comp.processes[i_proc],
            cb = tile.callbacks[name][i_cb],
            u = Strat.withaxes(u, tile),
            du = Strat.withaxes(get_du(integrator), tile),
            t = t;
            criterion(cb, layer, process, Strat.getstate(Val{name}(), tile, u, du, t))
        end
    end
end
function _affectfunc(::Val{name}, i_layer, i_proc, i_cb) where name
    function _affect!(integrator)
        let tile=Tile(integrator),
            comp = tile.strat[i_layer],
            layer = comp.layer,
            process = comp.processes[i_proc],
            cb = tile.callbacks[name][i_cb],
            u = Strat.withaxes(integrator.u, tile),
            du = Strat.withaxes(get_du(integrator), tile),
            t = integrator.t;
            affect!(cb, layer, process, Strat.getstate(Val{name}(), tile, u, du, t))
        end
    end
end
_diffeqcallback(::Discrete, ::Val{name}, i_layer, i_proc, i_cb) where name = DiffEqCallbacks.DiscreteCallback(
    _criterionfunc(Val{name}(), i_layer, i_proc, i_cb),
    _affectfunc(Val{name}(), i_layer, i_proc, i_cb),
    # todo: initialize and finalize?
)
_diffeqcallback(::Continuous, ::Val{name}, i_layer, i_proc, i_cb) where name = DiffEqCallbacks.ContinuousCallback(
    _criterionfunc(Val{name}(), i_layer, i_proc, i_cb),
    _affectfunc(Val{name}(), i_layer, i_proc, i_cb),
    # todo: initialize and finalize?
)
function _makecallbacks(component::StratComponent{L,P,name}, i_layer) where {L,P,name}
    cbs = []
    i_cb = 1
    for (i_proc, proc) in enumerate(component.processes)
        for callback in callbacks(component.layer, proc)
            push!(cbs, _diffeqcallback(CallbackStyle(callback), Val{name}(), i_layer, i_proc, i_cb))
            i_cb += 1
        end
    end
    return Tuple(cbs)
end
function _makecallbacks(tile::Tile)
    diffeq_callbacks = tuplejoin((_makecallbacks(comp, i) for (i,comp) in enumerate(tile.strat))...)
    return diffeq_callbacks
end

end
