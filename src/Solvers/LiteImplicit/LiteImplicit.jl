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

struct CGLiteSolution{TT,Tu<:AbstractVector{TT},Tt,Tprob} <: SciMLBase.AbstractODESolution{TT,1,Tu}
    prob::Tprob
    u::Vector{Tu}
    t::Vector{Tt}
end
Base.getproperty(sol::CGLiteSolution, name::Symbol) = name == :H ? getfield(sol, :u) : getfield(sol, name)
Base.propertynames(sol::CGLiteSolution) = (fieldnames(typoeof(sol))..., :H)
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

mutable struct CGLiteIntegrator{Talg,Tu,Tt,Tp,Tsol,Tcache} <: SciMLBase.DEIntegrator{Talg,true,Tu,Tt}
    alg::Talg
    cache::Tcache
    sol::Tsol
    u::Tu
    p::Tp
    t::Tt
    dt::Tt
    step::Int
end
SciMLBase.done(integrator::CGLiteIntegrator) = integrator.t >= integrator.sol.prob.tspan[end]
SciMLBase.get_du(integrator::CGLiteIntegrator) = integrator.cache.du

function DiffEqBase.__init(prob::CryoGridProblem, alg::LiteImplicitEuler, args...; dt=24*3600.0, kwargs...)
    tile = Tile(prob.f)
    grid = tile.grid
    u0 = copy(collect(prob.u0))
    nsteps = Int(ceil((prob.tspan[2] - prob.tspan[1]) / dt)) + 1
    u_storage = Vector{typeof(u0)}(undef, nsteps)
    t_storage = zeros(nsteps)
    # initialize storage
    u_storage[1] = u0
    t_storage[1] = prob.tspan[1]
    for i in 2:length(u_storage)
        u_storage[i] = similar(u0)
    end
    # reset SavedValues on tile.hist
    stateproto = prob.savefunc(tile, u0, similar(u0))
    savevals = SavedValues(Float64, typeof(stateproto))
    tile.hist.vals = savevals
    sol = CGLiteSolution(prob, u_storage, t_storage)
    cache = LiteImplicitEulerCache(
        similar(prob.u0), # should have ComponentArray type
        similar(prob.u0),
        similar(u0, length(prob.u0.H)),
        similar(u0, length(prob.u0.H)),
        similar(u0, length(prob.u0.H)),
        similar(u0, length(prob.u0.H)),
        similar(u0, length(prob.u0.H)),
        similar(u0, length(prob.u0.H)-1),
        similar(u0, length(prob.u0.H)),
        similar(u0, length(prob.u0.H)-1),
        similar(u0, length(prob.u0.H)),
    )
    return CGLiteIntegrator(alg, cache, sol, copy(prob.u0), collect(prob.p), prob.tspan[1], dt, 1)
end

function DiffEqBase.__solve(prob::CryoGridProblem, alg::LiteImplicitEuler, args...; dt=24*3600.0, kwargs...)
    integrator = __init(prob, alg, args...; dt=dt, kwargs...)
    for i in integrator end
    return integrator.sol
end

include("step.jl")

end
