module CryoGrid

using Base: @propagate_inbounds
using Interpolations
using TimeSeries
using Dates
using Lazy
using DataStructures: SortedDict, OrderedDict
using AxisArrays
using StaticArrays
using ComponentArrays
using DiffEqBase
using Reexport
@reexport using Unitful
@reexport using IntervalSets
@reexport using Parameters
@reexport using DifferentialEquations

# main version
include("core/core.jl")
include("layers/layers.jl")
include("processes/processes.jl")

# implicit version
# include("implicit/CryoGridImplicit.jl")

end # module
