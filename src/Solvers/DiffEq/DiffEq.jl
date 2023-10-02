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
using Requires
using Unitful
using UnPack

using DiffEqBase

# re-export DiffEqCallbacks
@reexport using DiffEqCallbacks

export TDMASolver
include("linsolve.jl")

export CryoGridEnsembleSetup, CryoGridEnsembleProblem
include("ensemble.jl")

function __init__()
    # OrdinaryDiffEq compatibility
    @require OrdinaryDiffEq="1dea7af3-3e70-54e6-95c3-0bf5283fa5ed" begin
        using .OrdinaryDiffEq
        using .OrdinaryDiffEq: NLSolver
        # re-export selected types from OrdinaryDiffEq;
        export OrdinaryDiffEq
        # explicit methods
        export Euler, Heun, DP5, Tsit5, ROCK2, ROCK4, SSPRK22, SSPRK33, SSPRK43
        # implicit methods
        export ImplicitEuler, ImplicitMidpoint, Trapezoid, SSPSDIRK2
        # nonlinear solvers
        export NLNewton, NLCGLite
        include("ode_solvers.jl")
    end
end

end
