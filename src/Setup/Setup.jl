module Setup

import ForwardDiff
import ReverseDiff
import Zygote

import CryoGrid: Layer, Top, Bottom, SubSurface, Process, SubSurfaceProcess, BoundaryProcess, System
import CryoGrid: variables, initialcondition!, prognosticstep!, diagnosticstep!, interact!, observe

using CryoGrid.InputOutput
using CryoGrid.Layers
using CryoGrid.Numerics
using CryoGrid.Processes
using CryoGrid.Utils

using ComponentArrays
using DataStructures: OrderedDict
using Dates
using DimensionalData
using IntervalSets
using Lazy: @>>, @>, groupby
using LinearAlgebra
using Unitful
using RecursiveArrayTools
using Reexport

@reexport using DiffEqBase: solve, init, ODEProblem, SciMLBase
@reexport using OrdinaryDiffEq
@reexport using SimulationLogs

export CryoGridSetup, CryoGridOutput, withaxes, getstate, getvar

include("stratigraphy.jl")
include("setup.jl")
include("problem.jl")
include("output.jl")

end