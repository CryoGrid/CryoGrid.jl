using DataStructures: SortedSet

abstract type CryoGridODEAlgorithm <: SciMLBase.AbstractODEAlgorithm end

abstract type AbstractCryoGridSolution{TT,N,Tu} <: SciMLBase.AbstractODESolution{TT,N,Tu} end

mutable struct CryoGridSolution{TT,Tu<:AbstractVector{TT},Tt,Talg,Tprob} <: AbstractCryoGridSolution{TT,1,Tu}
    prob::Tprob
    u::Vector{Tu}
    t::Vector{Tt}
    alg::Talg
    retcode::ReturnCode.T
end
Base.getproperty(sol::CryoGridSolution, name::Symbol) = name == :H ? getfield(sol, :u) : getfield(sol, name)
Base.propertynames(sol::CryoGridSolution) = (fieldnames(typeof(sol))..., :H)
function Base.show(io::IO, m::MIME"text/plain", sol::CryoGridSolution)
    println(io, string("retcode: ", sol.retcode))
    # TODO: inmplement interpolation interface?
    println(io, string("Interpolation: 1st order linear"))
    print(io, "t: ")
    show(io, m, sol.t)
    println(io)
    print(io, "u: ")
    show(io, m, sol.u)
end
# evaluate the solution at arbitrary time t
function (sol::CryoGridSolution)(t::Float64)
    N = length(sol.u[1])
    i = searchsortedfirst(sol.t, t)
    if i > length(sol.t)
        return sol.u[end]
    elseif i == 1
        return sol.u[1]
    elseif t == sol.t[i]
        return sol.u[i]
    else
        # simple linear initerpolation between tᵢ₋₁ and tᵢ
        t₁ = sol.t[i-1]
        t₂ = sol.t[i]
        u₁ = sol.u[i-1]
        u₂ = sol.u[i]
        return sol.u[i-1] .+ t.*(u₂ - u₁) ./ (t₂ - t₁)
    end
end
function DiffEqBase.sensitivity_solution(sol::CryoGridSolution, u, t)
    T = eltype(eltype(u))
    N = length((size(sol.prob.u0)..., length(u)))
    interp = LinearInterpolation(t, u)
    ODESolution{T, N}(u, nothing, nothing, t,
                      nothing, sol.prob,
                      sol.alg, interp,
                      true, length(sol.t),
                      nothing, sol.retcode)
end

Base.@kwdef mutable struct CryoGridIntegratorOptions{Tt,Tsaveat,Ttstops}
    dtmax::Tt = 3600.0
    dtmin::Tt = one(typeof(dtmax))
    saveat::Tsaveat = typeof(dtmax)[]
    tstops::Ttstops = SortedSet{typeof(dtmin)}()
    save_everystep::Bool = true
    just_hit_tstop::Bool = false
    stop_at_next_tstop::Bool = true
end

mutable struct CryoGridIntegrator{Talg,Tu,Tt,Tp,Topts,Tsol,Tcache} <: SciMLBase.AbstractODEIntegrator{Talg,true,Tu,Tt}
    alg::Talg
    cache::Tcache
    opts::Topts
    sol::Tsol
    u::Tu
    p::Tp
    t::Tt
    dt::Tt
    tdir::Int
    step::Int
end
SciMLBase.done(integrator::CryoGridIntegrator) = integrator.t >= integrator.sol.prob.tspan[end]
SciMLBase.get_du(integrator::CryoGridIntegrator) = integrator.cache.du
SciMLBase.add_tstop!(integrator::CryoGridIntegrator, t) = push!(integrator.opts.tstops, t)
SciMLBase.postamble!(integrator::CryoGridIntegrator) = nothing

# add tstop by default because we don't support fancy interpolation
DiffEqBase.step!(integrator::CryoGridIntegrator, dt) = step!(integrator, dt, true)

function CommonSolve.step!(integrator::CryoGridIntegrator)
    handle_tstops!(integrator)
    perform_step!(integrator)
    saveat!(integrator)    
    # ajust dt based on options
    integrator.dt = min(integrator.opts.dtmax, integrator.dt)
end

function CommonSolve.solve!(integrator::CryoGridIntegrator)
    for i in integrator end
    if integrator.sol.retcode == ReturnCode.Default
        integrator.sol.retcode = ReturnCode.Success
    end
    # if no save points are specified, save final state
    if isempty(integrator.sol.prob.saveat)
        prob.savefunc(integrator.u, integrator.t, integrator)
        push!(integrator.sol.u, integrator.u)
        push!(integrator.sol.t, integrator.t)
    end
    return integrator.sol
end

# CryoGridIntegrator interface

perform_step!(integrator::CryoGridIntegrator) = error("perform_step! not implemented for algorithm $(integrator.alg)")

function saveat!(integrator::CryoGridIntegrator)
    prob = integrator.sol.prob
    saveat = integrator.opts.saveat
    t_saves = integrator.sol.t
    u_saves = integrator.sol.u
    res = searchsorted(saveat, integrator.t)
    i_next = first(res)
    i_prev = last(res)
    dtsave = if i_next == i_prev
        prob.savefunc(integrator.u, integrator.t, integrator)
        push!(u_saves, copy(integrator.u))
        push!(t_saves, integrator.t)
        Inf
    elseif i_next > length(saveat)
        Inf
    else
        saveat[i_next] - integrator.t
    end
    integrator.dt = min(integrator.dt, dtsave)
end

function handle_tstops!(integrator::CryoGridIntegrator)
    if !isempty(integrator.opts.tstops)
        next_tstop = first(integrator.opts.tstops)
        dt_to_stop = next_tstop - integrator.t
        if dt_to_stop > zero(dt_to_stop)
            integrator.dt = min(integrator.dt, dt_to_stop)
        else
            pop!(integrator.opts.tstops)
        end
    end
end

expandtstep(tstep::Number, tspan) = tspan[1]:tstep:tspan[end]
expandtstep(tstep::AbstractVector, tspan) = tstep

# Output

"""
    CryoGridOutput(sol::AbstractCryoGridSolution, tspan::NTuple{2,Float64}=(-Inf,Inf))

Constructs a `CryoGridOutput` from the given `ODESolution`. Optional argument `tspan` restricts the time span of the output.
"""
InputOutput.CryoGridOutput(sol::AbstractCryoGridSolution, tspan::NTuple{2,DateTime}) = CryoGridOutput(sol, convert_tspan(tspan))
function InputOutput.CryoGridOutput(sol::AbstractCryoGridSolution, tspan::NTuple{2,Float64}=(-Inf,Inf))
    # Helper functions for mapping variables to appropriate DimArrays by grid/shape.
    withdims(var::Var{name,<:CryoGrid.OnGrid{Cells}}, arr, grid, ts) where {name} = DimArray(arr*one(vartype(var))*varunits(var), (Z(round.(typeof(1.0u"m"), cells(grid), digits=5)),Ti(ts)))
    withdims(var::Var{name,<:CryoGrid.OnGrid{Edges}}, arr, grid, ts) where {name} = DimArray(arr*one(vartype(var))*varunits(var), (Z(round.(typeof(1.0u"m"), edges(grid), digits=5)),Ti(ts)))
    withdims(var::Var{name}, arr, zs, ts) where {name} = DimArray(arr*one(vartype(var))*varunits(var), (Dim{name}(1:size(arr,1)),Ti(ts)))
    save_interval = ClosedInterval(tspan...)
    tile = Tile(sol.prob.f) # Tile
    grid = Grid(tile.grid.*u"m")
    savecache = sol.prob.savefunc.cache
    ts = savecache.t # use save cache time points
    # check if last value is duplicated
    ts = ts[end] == ts[end-1] ? ts[1:end-1] : ts
    t_mask = map(∈(save_interval), ts) # indices within t interval
    u_all = sol.(ts[t_mask])
    u_mat = reduce(hcat, u_all) # build prognostic state from continuous solution
    pax = ComponentArrays.indexmap(getaxes(tile.state.uproto)[1])
    # get saved diagnostic states and timestamps only in given interval
    savedstates = savecache.vals[1:length(ts)][t_mask]
    ts_datetime = Dates.epochms2datetime.(round.(ts[t_mask]*1000.0))
    allvars = variables(tile)
    progvars = tuplejoin(filter(isprognostic, allvars), filter(isalgebraic, allvars))
    diagvars = filter(isdiagnostic, allvars)
    fluxvars = filter(isflux, allvars)
    outputs = OrderedDict()
    # add all on-grid prognostic variables
    for var in filter(isongrid, progvars)
        name = varname(var)
        outputs[name] = withdims(var, u_mat[pax[name],:], grid, ts_datetime)
    end
    # add all on-grid diagnostic variables
    for var in filter(isongrid, tuplejoin(diagvars, fluxvars))
        name = varname(var)
        states = collect(skipmissing([name ∈ keys(state) ? state[name] : missing for state in savedstates]))
        if length(states) == length(ts_datetime)
            arr = reduce(hcat, states)
            outputs[name] = withdims(var, arr, grid, ts_datetime)
        end
    end
    # handle per-layer variables
    for layer in layernames(tile.strat)
        # if layer name appears in saved states or prognostic state axes, then add these variables to the output.
        if haskey(savedstates[1], layer) || haskey(pax, layer)
            # map over all savedstates and create named tuples for each time step
            layerouts = map(u_all, savedstates) do u, state
                layerout = OrderedDict()
                if haskey(state, layer)
                    layerstate = state[layer]
                    for var in keys(layerstate)
                        layerout[var] = layerstate[var]
                    end
                else
                end
                # convert to named tuple
                diagnostic_output = (;layerout...)
                if haskey(u, layer)
                    u_layer = u[layer]
                    prognostic_output = (;map(name -> name => u_layer[name], keys(u_layer))...)
                    return merge(prognostic_output, diagnostic_output)
                else
                    return diagnostic_output
                end
            end
            layerouts_combined = reduce(layerouts[2:end]; init=layerouts[1]) do out1, out2
                map(vcat, out1, out2)
            end
            # for each variable in the named tuple, find the corresponding variables
            layervars = (; map(name -> name => first(filter(var -> varname(var) == name, allvars)), keys(layerouts_combined))...)
            # map each output to a variable and call withdims to wrap in a DimArray
            outputs[layer] = map((var,out) -> withdims(var, reshape(out,1,:), nothing, ts_datetime), layervars, layerouts_combined)
        end
    end
    return CryoGridOutput(ts_datetime, sol, (;outputs...))
end
