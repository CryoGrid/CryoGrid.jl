"""
Driver module SciML diffeq solvers.
"""
module DiffEq

using CryoGrid
using CryoGrid: Event, ContinuousEvent, DiscreteEvent, ContinuousTrigger, Increasing, Decreasing
using CryoGrid.Drivers
using CryoGrid.InputOutput
using CryoGrid.Numerics
using CryoGrid.Physics: Heat
using CryoGrid.Strat: Stratigraphy, StratComponent
using CryoGrid.Utils

import CryoGrid.Strat

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

using DiffEqBase
using DiffEqBase.SciMLBase
using DiffEqCallbacks

import DiffEqCallbacks

@reexport using OrdinaryDiffEq
@reexport using DiffEqBase: solve, init, ODEProblem, SciMLBase

export CFLStepLimiter
include("steplimiters.jl")
export TDMASolver
include("solvers.jl")
include("callbacks.jl")
export CryoGridProblem
include("problem.jl")
include("output.jl")

# Add method dispatches for other CryoGrid methods on DEIntegrator type
"""
    Tile(integrator::SciMLBase.DEIntegrator)

Constructs a `Tile` from a `SciMLBase` integrator.
"""
function Strat.Tile(integrator::SciMLBase.DEIntegrator)
    tile = integrator.sol.prob.f.f
    return Strat.updateparams(tile, Strat.withaxes(integrator.u, tile), integrator.p, integrator.t)
end
"""
    getstate(integrator::SciMLBase.DEIntegrator)
    getstate(layername::Symbol, integrator::SciMLBase.DEIntegrator)

Builds the state named tuple for `layername` given an initialized integrator.
"""
Strat.getstate(integrator::SciMLBase.DEIntegrator) = Strat.getstate(Tile(integrator), integrator.u, get_du(integrator), integrator.t)
Strat.getstate(layername::Symbol, integrator::SciMLBase.DEIntegrator) = Strat.getstate(Val{layername}(), integrator)
Strat.getstate(::Val{layername}, integrator::SciMLBase.DEIntegrator) where {layername} = Strat.getstate(Val{layername}(), Tile(integrator), integrator.u, get_du(integrator), integrator.t)
"""
    getvar(var::Symbol, integrator::SciMLBase.DEIntegrator)
"""
Numerics.getvar(var::Symbol, integrator::SciMLBase.DEIntegrator) = Numerics.getvar(Val{var}(), Tile(integrator), integrator.u)

end
