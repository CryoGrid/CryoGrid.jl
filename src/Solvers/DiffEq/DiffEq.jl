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
using SciMLBase
using DiffEqCallbacks

@reexport using OrdinaryDiffEq
@reexport using DiffEqBase: solve, init, ODEProblem

export TDMASolver
include("linsolve.jl")

export NLCGLite
include("nlsolve/nlsolve.jl")

export CryoGridEnsembleSetup, CryoGridEnsembleProblem
include("ensemble.jl")

# solve/init interface
function DiffEqBase.__solve(prob::CryoGridProblem, alg::Union{OrdinaryDiffEqAlgorithm, OrdinaryDiffEq.DAEAlgorithm}, args...; kwargs...)
    ode_prob = ODEProblem(prob)
    return DiffEqBase.solve(ode_prob, alg, args...; kwargs...)
end
function DiffEqBase.__init(prob::CryoGridProblem, alg::Union{OrdinaryDiffEqAlgorithm, OrdinaryDiffEq.DAEAlgorithm}, args...; kwargs...)
    ode_prob = ODEProblem(prob)
    return DiffEqBase.init(ode_prob, alg, args...; kwargs...)
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
function Tiles.Tile(f::ODEFunction)
    extract_f(tile::Tile) = tile
    extract_f(f::ODEFunction) = f.f
    extract_f(f::DiffEqBase.Void) = f.f
    extract_f(f) = SciMLBase.unwrapped_f(f)
    return extract_f(f.f)
end
"""
    getstate(integrator::SciMLBase.DEIntegrator)
    getstate(layername::Symbol, integrator::SciMLBase.DEIntegrator)

Builds the state named tuple for `layername` given an initialized integrator.
"""
Tiles.getstate(integrator::SciMLBase.DEIntegrator) = Tiles.getstate(Tile(integrator), integrator.u, get_du(integrator), integrator.t)
Tiles.getstate(layername::Symbol, integrator::SciMLBase.DEIntegrator) = Tiles.getstate(Val{layername}(), integrator)
Tiles.getstate(::Val{layername}, integrator::SciMLBase.DEIntegrator) where {layername} = Tiles.getstate(Val{layername}(), Tile(integrator), integrator.u, get_du(integrator), integrator.t)
"""
    getvar(var::Symbol, integrator::SciMLBase.DEIntegrator)
"""
Numerics.getvar(var::Symbol, integrator::SciMLBase.DEIntegrator) = Numerics.getvar(Val{var}(), Tile(integrator), integrator.u)

end
