module LiteImplicit

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils

using DiffEqBase, DiffEqCallbacks
using Interpolations

import DiffEqBase
import SciMLBase

export LiteImplicitEuler

Base.@kwdef struct LiteImplicitEuler <: SciMLBase.AbstractODEAlgorithm
    miniters::Int = 2
    maxiters::Int = 1000
    tolerance::Float64 = 1e-3
    verbose::Bool = true
end
DiffEqBase.check_prob_alg_pairing(::CryoGridProblem, ::LiteImplicitEuler) = nothing
DiffEqBase.check_prob_alg_pairing(prob, alg::LiteImplicitEuler) = throw(DiffEqBase.ProblemSolverPairingError(prob, alg))

struct LiteImplicitEulerCache{Tu,TA} <: SciMLBase.DECache
    uprev::Tu
    du::Tu
    H::TA
    T_new::TA
    resid::TA
    Sc::TA
    Sp::TA
    A::TA
    B::TA
    C::TA
    D::TA
end

mutable struct CGLiteSolution{TT,Tu<:AbstractVector{TT},Tt,Talg,Tprob} <: SciMLBase.AbstractODESolution{TT,1,Tu}
    prob::Tprob
    u::Vector{Tu}
    t::Vector{Tt}
    alg::Talg
    retcode::Symbol
end
Base.getproperty(sol::CGLiteSolution, name::Symbol) = name == :H ? getfield(sol, :u) : getfield(sol, name)
Base.propertynames(sol::CGLiteSolution) = (fieldnames(typeof(sol))..., :H)
function Base.show(io::IO, m::MIME"text/plain", sol::CGLiteSolution)
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
function (sol::CGLiteSolution)(t::Float64)
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
function DiffEqBase.sensitivity_solution(sol::CGLiteSolution, u, t)
    T = eltype(eltype(u))
    N = length((size(sol.prob.u0)..., length(u)))
    interp = LinearInterpolation(t, u)
    ODESolution{T, N}(u, nothing, nothing, t,
                      nothing, sol.prob,
                      sol.alg, interp,
                      true, length(sol.t),
                      nothing, sol.retcode)
end

mutable struct CGLiteIntegrator{Talg,Tu,Tt,Tp,Tsol,Tcache} <: SciMLBase.DEIntegrator{Talg,true,Tu,Tt}
    alg::Talg
    cache::Tcache
    sol::Tsol
    u::Tu
    p::Tp
    t::Tt
    dt::Tt
    tdir::Int
    step::Int
end
SciMLBase.done(integrator::CGLiteIntegrator) = integrator.t >= integrator.sol.prob.tspan[end]
SciMLBase.get_du(integrator::CGLiteIntegrator) = integrator.cache.du

function DiffEqBase.__init(prob::CryoGridProblem, alg::LiteImplicitEuler, args...; dt=24*3600.0, kwargs...)
    tile = Tile(prob.f)
    grid = tile.grid
    u0 = copy(prob.u0)
    t0 = prob.tspan[1]
    nsteps = Int(ceil((prob.tspan[2] - t0) / dt)) + 1
    # initialize storage
    u_storage = [u0]
    t_storage = [prob.tspan[1]]
    # evaluate tile at initial condition
    tile = Strat.resolve(Tile(prob.f), u0, prob.p, t0)
    CryoGrid.Strat.step!(tile, zero(u0), u0, prob.p, t0, dt)
    # reset SavedValues on tile.hist
    initialsave = prob.savefunc(tile, u0, similar(u0))
    savevals = SavedValues(Float64, typeof(initialsave))
    push!(savevals.saveval, initialsave)
    push!(savevals.t, t0)
    tile.hist.vals = savevals
    sol = CGLiteSolution(prob, u_storage, t_storage, alg, :Default)
    cache = LiteImplicitEulerCache(
        similar(prob.u0), # should have ComponentArray type
        similar(prob.u0),
        similar(u0, eltype(u0), length(prob.u0.H)),
        similar(u0, eltype(u0), length(prob.u0.H)),
        similar(u0, eltype(u0), length(prob.u0.H)),
        similar(u0, eltype(u0), length(prob.u0.H)),
        similar(u0, eltype(u0), length(prob.u0.H)),
        similar(u0, eltype(u0), length(prob.u0.H)-1),
        similar(u0, eltype(u0), length(prob.u0.H)),
        similar(u0, eltype(u0), length(prob.u0.H)-1),
        similar(u0, eltype(u0), length(prob.u0.H)),
    )
    p = isnothing(prob.p) ? prob.p : collect(prob.p)
    return CGLiteIntegrator(alg, cache, sol, u0, p, t0, convert(eltype(prob.tspan), dt), 1, 1)
end

function DiffEqBase.__solve(prob::CryoGridProblem, alg::LiteImplicitEuler, args...; dt=24*3600.0, kwargs...)
    integrator = DiffEqBase.__init(prob, alg, args...; dt=dt, kwargs...)
    for i in integrator end
    if integrator.sol.retcode == :Default
        integrator.sol.retcode = :Success
    end
    return integrator.sol
end

include("step.jl")

end
