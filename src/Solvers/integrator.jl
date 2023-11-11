using DataStructures: SortedSet

abstract type CryoGridODEAlgorithm <: SciMLBase.AbstractODEAlgorithm end

mutable struct CryoGridSolution{TT,Tu<:AbstractVector{TT},Tt,Talg,Tprob} <: SciMLBase.AbstractODESolution{TT,1,Tu}
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
    dtmax::Tt = 24*3600.0
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

function DiffEqBase.step!(integrator::CryoGridIntegrator)
    handle_tstops!(integrator)
    perform_step!(integrator)
    saveat!(integrator)    
    # ajust dt based on options
    integrator.dt = min(integrator.opts.dtmax, integrator.dt)
end

function DiffEqBase.__solve(prob::CryoGridProblem, alg::CryoGridODEAlgorithm, args...; kwargs...)
    integrator = DiffEqBase.__init(prob, alg, args...; kwargs...)
    for i in integrator end
    if integrator.sol.retcode == ReturnCode.Default
        integrator.sol.retcode = ReturnCode.Success
    end
    # if no save points are specified, save final state
    if isempty(prob.saveat)
        tile = Tile(integrator)
        push!(tile.data.outputs.saveval, prob.savefunc(tile, integrator.u, get_du(integrator)))
        push!(tile.data.outputs.t, ForwardDiff.value(integrator.t))
        push!(integrator.sol.u, integrator.u)
        push!(integrator.sol.t, integrator.t)
    end
    return integrator.sol
end

perform_step!(integrator::CryoGridIntegrator) = error("perform_step! not implemented for algorithm $(integrator.alg)")

function saveat!(integrator::CryoGridIntegrator)
    tile = Tile(integrator)
    du = get_du(integrator)
    prob = integrator.sol.prob
    saveat = prob.saveat
    t_saves = integrator.sol.t
    u_saves = integrator.sol.u
    res = searchsorted(saveat, integrator.t)
    i_next = first(res)
    i_prev = last(res)
    dtsave = if i_next == i_prev    
        push!(tile.data.outputs.saveval, integrator.sol.prob.savefunc(tile, integrator.u, du))
        push!(tile.data.outputs.t, ForwardDiff.value(integrator.t))
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
