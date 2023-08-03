"""
Driver module SciML diffeq solvers.
"""
module DiffEq

using CryoGrid
using CryoGrid: Event, ContinuousEvent, DiscreteEvent, ContinuousTrigger, Increasing, Decreasing
using CryoGrid.InputOutput
using CryoGrid.Numerics
using CryoGrid.Utils

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
using UnPack

using DiffEqBase
using OrdinaryDiffEq

# re-export DiffEqCallbacks
@reexport using DiffEqCallbacks

# re-export selected types from OrdinaryDiffEq;
export OrdinaryDiffEq
# explicit methods
export Euler, Heun, DP5, Tsit5, ROCK2, ROCK4, SSPRK22, SSPRK33, SSPRK43
# implicit methods
export ImplicitEuler, ImplicitMidpoint, Trapezoid, SSPSDIRK2

export TDMASolver
include("linsolve.jl")

export NLCGLite
include("nlsolve/nlsolve.jl")

export CryoGridEnsembleSetup, CryoGridEnsembleProblem
include("ensemble.jl")

# solve/init interface
function DiffEqBase.__solve(prob::CryoGridProblem, alg::Union{OrdinaryDiffEqAlgorithm, OrdinaryDiffEq.DAEAlgorithm}, args...; saveat=prob.saveat, kwargs...)
    ode_prob = ODEProblem(prob)
    return DiffEqBase.solve(ode_prob, alg, args...; saveat, kwargs...)
end
function DiffEqBase.__init(prob::CryoGridProblem, alg::Union{OrdinaryDiffEqAlgorithm, OrdinaryDiffEq.DAEAlgorithm}, args...; saveat=prob.saveat, kwargs...)
    ode_prob = ODEProblem(prob)
    return DiffEqBase.init(ode_prob, alg, args...; saveat, kwargs...)
end

# Add method dispatches for other CryoGrid methods on DEIntegrator type
"""
    Tile(integrator::SciMLBase.DEIntegrator)

Constructs a `Tile` from a `SciMLBase` integrator.
"""
function Tiles.Tile(integrator::SciMLBase.DEIntegrator)
    tile = Tiles.Tile(integrator.sol.prob.f)
    du = get_du(integrator)
    u = integrator.u
    return Tiles.resolve(tile, Tiles.withaxes(u, tile), integrator.p, integrator.t)
end

end
