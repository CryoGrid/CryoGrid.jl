"""
    Common

Common types and functions used by all other CryoGrid.jl modules.
"""
module Common

import Base.==

using ArraysOfArrays
using Base: @propagate_inbounds, Float64
using ComponentArrays
using DataStructures: SortedDict
using Interpolations: Interpolations, Gridded, Linear, Flat, Line, interpolate, extrapolate
using Lazy: @>>, @>
using Reexport
using TimeSeries

# Re-exported packages
@reexport using Dates
@reexport using DiffEqBase
@reexport using DiffEqCallbacks
@reexport using DimensionalData: DimArray, X, Y, Z, Dim, dims
@reexport using IntervalSets
@reexport using OrdinaryDiffEq
@reexport using Parameters
@reexport using Unitful
@reexport using SimulationLogs

# Submodules
include("utils/Utils.jl")
include("math/Math.jl")

@reexport using .Utils
@reexport using .Math

# Common components
include("types.jl")
include("grid.jl")
include("forcing.jl")
include("variables.jl")
include("stratigraphy.jl")

end