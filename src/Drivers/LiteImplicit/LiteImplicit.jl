module LiteImplicit

using CryoGrid.Numerics
using CryoGrid.Utils
using CryoGrid.Strat

using DiffEqBase

import DiffEqBase: __solve, __init, step!
import SciMLBase: done

struct CryoGridLiteProblem{Tu,Tt,Tp,TT} <: SciMLBase.AbstractODEProblem{Tu,Tt,true}
    tile::TT
    u0::Tu
    p::Tp
    tspan::NTuple{2,Tt}
end

Base.@kwdef struct LiteImplicitEuler <: SciMLBase.AbstractDEAlgorithm
    miniters::Int = 10
    maxiters::Int = 1000
    tolerance::Float64 = 1e-3
end

struct LiteImplicitCache{TA} <: SciMLBase.DECache
    H::TA
    dH::TA
    T_new::TA
    resid::TA
    dxn::TA
    Sc::TA
    Sp::TA
    bp::TA
    ap::TA
    ans::TA
    A::TA
    B::TA
    C::TA
    D::TA
end

struct CGLiteSolution{T,Tu,Tt,Tprob} <: SciMLBase.AbstractODESolution{T,1,Tu}
    u::Vector{Tu}
    t::Vector{Tt}
    prob::Tprob
end

mutable struct CGLiteIntegrator{Talg,Tu,Tt,Tp,Tsol,Tcache} <: SciMLBase.DEIntegrator{Talg,true,Tu,Tt}
    alg::Talg
    cache::Tcache
    sol::Tsol
    u::Tu
    p::Tp
    t::Tt
    dt::Tt
end
done(integrator::CGLiteIntegrator) = integrator.t >= integrator.sol.prob.tspan[end]

function __init(prob::CryoGridLiteProblem, alg::LiteImplicitEuler, args...; dt=24*3600.0, kwargs...)
    sol = CGLiteSolution([copy(prob.u0)], [prob.tspan[1]], prob)
    grid = prob.tile.grid
    cache = LiteImplicitCache(
        copy(prob.u0),
        zero(prob.u0),
        zero(prob.u0),
        zero(prob.u0),
        zero(prob.u0),
        zero(prob.u0),
        similar(prob.u0, length(cells(grid))-1),
        zero(prob.u0),
        zero(prob.u0),
        zero(prob.u0),
        similar(prob.u0, length(prob.u0)-1),
        similar(prob.u0),
        similar(prob.u0, length(prob.u0)-1),
        similar(prob.u0)
    )
    return CGLiteIntegrator(alg, cache, sol, copy(prob.u0), prob.p, prob.tspan[1], dt)
end

function __solve(prob::CryoGridLiteProblem, alg::LiteImplicitEuler, args...; dt=24*3600.0, kwargs...)
    integrator = __init(prob, alg, args...; dt=24*3600.0, kwargs...)
    for i in integrator end
    return integrator.sol
end

end
