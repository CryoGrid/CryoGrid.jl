"""
    CryoGridDiffEqSolution

Wrapper types for ODE solutions from `OrdinaryDiffEq`.
"""
struct CryoGridDiffEqSolution{
    TT,
    N,
    Tu,
    solType<:SciMLBase.AbstractODESolution{TT,N,Tu}
} <: CryoGrid.AbstractCryoGridSolution{TT,N,Tu}
    prob::CryoGridProblem
    sol::solType
end

Base.propertynames(sol::CryoGridDiffEqSolution) = (
    :prob,
    Base.propertynames(sol.sol)...
)

function Base.getproperty(sol::CryoGridDiffEqSolution, name::Symbol)
    if name == :prob
        return getfield(sol, name)
    else
        return getproperty(getfield(sol, :sol), name)
    end
end

"""
    CryoGridDiffEqIntegrator

Wrapper type for ODE integrators from `OrdinaryDiffEq`.
"""
struct CryoGridDiffEqIntegrator{
    algType,
    Tu,
    Tt,
    intType<:SciMLBase.AbstractODEIntegrator{algType,true,Tu,Tt}
} <: SciMLBase.AbstractODEIntegrator{algType,true,Tu,Tt}
    prob::CryoGridProblem
    integrator::intType
end

SciMLBase.done(solver::CryoGridDiffEqIntegrator) = SciMLBase.done(solver.integrator)
SciMLBase.get_du(integrator::CryoGridDiffEqIntegrator) = get_du(integrator.integrator)
SciMLBase.add_tstop!(integrator::CryoGridDiffEqIntegrator, t) = push!(integrator.opts.tstops, t)
SciMLBase.postamble!(solver::CryoGridDiffEqIntegrator) = SciMLBase.postamble!(solver.integrator)

function Base.show(io::IO, mime::MIME"text/plain", integrator::CryoGridDiffEqIntegrator)
    println(io, "CryoGrid ODE integrator:")
    show(io, mime, integrator.integrator)
end

Base.propertynames(integrator::CryoGridDiffEqIntegrator) = (
    :prob,
    :integrator,
    Base.propertynames(integrator.integrator)...
)

function Base.getproperty(integrator::CryoGridDiffEqIntegrator, name::Symbol)
    inner = getfield(integrator, :integrator)
    if name ∈ (:prob, :integrator)
        return getfield(integrator, name)
    elseif name == :sol
        return CryoGridDiffEqSolution(integrator.prob, inner.sol)
    else
        return getproperty(inner, name)
    end
end

function Base.setproperty!(integrator::CryoGridDiffEqIntegrator, name::Symbol, value)
    return setproperty!(getfield(integrator, :integrator), name, value)
end

# CommonSolve solve/init interface

function CommonSolve.init(
    prob::CryoGridProblem,
    alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm},
    args...;
    kwargs...
)
    ode_prob = ODEProblem(prob)
    ode_integrator = init(ode_prob, alg, args...; kwargs...)
    integrator = CryoGridDiffEqIntegrator(prob, ode_integrator)
    return integrator
end

function CommonSolve.step!(integrator::CryoGridDiffEqIntegrator, args...; kwargs...)
    rv = step!(integrator.integrator, args...; kwargs...)
    return rv
end

function CommonSolve.solve!(integrator::CryoGridDiffEqIntegrator)
    ode_sol = solve!(integrator.integrator)
    return CryoGridDiffEqSolution(integrator.prob, ode_sol)
end

# custom nonlinear solvers
include("nlsolve/nlsolve.jl")
