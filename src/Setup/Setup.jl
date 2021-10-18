module Setup

import ForwardDiff
import ReverseDiff

import CryoGrid: Layer, Top, Bottom, SubSurface, Process, SubSurfaceProcess, BoundaryProcess, CompoundProcess
import CryoGrid: variables, initialcondition!, prognosticstep!, diagnosticstep!, interact!, observe

using CryoGrid.InputOutput
using CryoGrid.Layers
using CryoGrid.Numerics
using CryoGrid.Processes
using CryoGrid.Utils

using ComponentArrays
using ConstructionBase
using DataStructures: OrderedDict
using Dates
using DiffEqCallbacks
using DimensionalData
using IfElse
using IntervalSets
using Lazy: @>>, @>, groupby
using LinearAlgebra
using Parameters
using PreallocationTools
using Unitful
using RecursiveArrayTools
using Reexport
using Setfield

import Flatten

@reexport using DiffEqBase: solve, init, ODEProblem, SciMLBase
@reexport using ModelParameters
@reexport using OrdinaryDiffEq
@reexport using SimulationLogs

export Stratigraphy, @Stratigraphy
export top, bottom, subsurface
export StratComponent, componentname, copmonenttypes, components, boundaries

export CryoGridSetup, CryoGridOutput, withaxes, getstate, getvar
export ParameterVector, LinearTrend, parameterize

include("stratigraphy.jl")
include("setup.jl")
include("parameterize.jl")
include("problem.jl")
include("output.jl")

end