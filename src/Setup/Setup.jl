module Setup

import ForwardDiff
import ReverseDiff
import Zygote

import CryoGrid.Interface: Top, Bottom

using CryoGrid.Interface
using CryoGrid.Layers
using CryoGrid.Numerics
using CryoGrid.Processes
using CryoGrid.Utils

using ArraysOfArrays
using ComponentArrays
using DataStructures: OrderedDict
using Dates
using DiffEqBase: ODEProblem, SciMLBase
using DimensionalData
using IntervalSets
using Lazy: @>>, @>, groupby
using LinearAlgebra
using Unitful
using RecursiveArrayTools
using Reexport

@reexport using OrdinaryDiffEq
@reexport using SimulationLogs

include("stratigraphy.jl")
include("setup.jl")
include("problem.jl")
include("output.jl")

end