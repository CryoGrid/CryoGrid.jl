@reexport using SimulationLogs

# solve/init interface
function DiffEqBase.__solve(prob::CryoGridProblem, alg::Union{OrdinaryDiffEqAlgorithm, OrdinaryDiffEq.DAEAlgorithm}, args...; saveat=prob.saveat, kwargs...)
    ode_prob = ODEProblem(prob)
    return DiffEqBase.solve(ode_prob, alg, args...; saveat, kwargs...)
end
function DiffEqBase.__init(prob::CryoGridProblem, alg::Union{OrdinaryDiffEqAlgorithm, OrdinaryDiffEq.DAEAlgorithm}, args...; saveat=prob.saveat, kwargs...)
    ode_prob = ODEProblem(prob)
    return DiffEqBase.init(ode_prob, alg, args...; saveat, kwargs...)
end

# custom nonlinear solvers
include("nlsolve/nlsolve.jl")
