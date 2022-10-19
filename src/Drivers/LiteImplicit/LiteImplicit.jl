module LiteImplicit

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils

using DiffEqBase
using SciMLBase

import DiffEqBase: __solve, __init, step!
import SciMLBase: done

export CryoGridLiteProblem, LiteImplicitEuler

struct CryoGridLiteProblem{Tu,Tt,Tp,TT,Tkw} <: SciMLBase.AbstractODEProblem{Tu,Tt,true}
    f::TT
    u0::Tu
    p::Tp
    tspan::NTuple{2,Tt}
    kwargs::Tkw
    function CryoGridLiteProblem(tile::TT, u0::Tu, tspan::NTuple{2}, p::Tp=p; kwargs...) where {TT<:HeatOnlyTile,Tu,Tp}
        tspan = convert_tspan(tspan)
        f = ODEFunction(tile)
        return new{Tu,eltype(tspan),Tp,typeof(f),typeof(kwargs)}(f, u0, p, tspan, kwargs)
    end
end

Base.@kwdef struct LiteImplicitEuler <: SciMLBase.AbstractDEAlgorithm
    miniters::Int = 2
    maxiters::Int = 1000
    tolerance::Float64 = 1e-3
end

struct LiteImplicitEulerCache{TA} <: SciMLBase.DECache
    H::TA
    dH::TA
    T_new::TA
    resid::TA
    dx::TA
    Sc::TA
    Sp::TA
    bp::TA
    ap::TA
    an::TA
    as::TA
    A::TA
    B::TA
    C::TA
    D::TA
end

struct CGLiteSolution{TT,Tu<:AbstractVector{TT},Tt,Tprob} <: SciMLBase.AbstractODESolution{TT,1,Tu}
    prob::Tprob
    u::Vector{Tu}
    t::Vector{Tt}
    T::Vector{Tu}
end
Base.getproperty(sol::CGLiteSolution, name::Symbol) = name == :H ? getfield(sol, :u) : getfield(sol, name)
Base.propertynames(sol::CGLiteSolution) = (fieldnames(typoeof(sol))..., :H)

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
done(integrator::CGLiteIntegrator) = integrator.t >= integrator.sol.prob.tspan[end]

function __init(prob::CryoGridLiteProblem, alg::LiteImplicitEuler, args...; dt=24*3600.0, kwargs...)
    tile = Tile(prob.f)
    grid = tile.grid
    u0 = copy(collect(prob.u0))
    T0 = similar(u0)
    # evaluate the step at the initial state to ensure correct starting values
    tile(similar(u0), u0, prob.p, prob.tspan[1])
    copyto!(T0, getvar(:T, tile, u0; interp=false))
    nsteps = Int(ceil((prob.tspan[2] - prob.tspan[1]) / dt)) + 1
    u_storage = Vector{typeof(u0)}(undef, nsteps)
    T_storage = Vector{typeof(T0)}(undef, nsteps)
    t_storage = zeros(nsteps)
    # initialize storage
    u_storage[1] = u0
    T_storage[1] = T0
    t_storage[1] = prob.tspan[1]
    for i in 2:length(u_storage)
        u_storage[i] = similar(u0)
        T_storage[i] = similar(T0)
    end
    sol = CGLiteSolution(prob, u_storage, t_storage, T_storage)
    cache = LiteImplicitEulerCache(
        copy(u0),
        zero(u0),
        zero(u0),
        zero(u0),
        similar(u0, length(cells(grid))-1),
        zero(u0),
        zero(u0),
        zero(u0),
        zero(u0),
        similar(u0, length(u0)-1),
        similar(u0, length(u0)-1),
        similar(u0, length(u0)-1),
        similar(u0),
        similar(u0, length(u0)-1),
        similar(u0)
    )
    return CGLiteIntegrator(alg, cache, sol, copy(prob.u0), collect(prob.p), prob.tspan[1], dt, 1)
end

function __solve(prob::CryoGridLiteProblem, alg::LiteImplicitEuler, args...; dt=24*3600.0, kwargs...)
    integrator = __init(prob, alg, args...; dt=dt, kwargs...)
    for i in integrator end
    return integrator.sol
end

include("step.jl")

end
