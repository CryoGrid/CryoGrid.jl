module CryoGrid

using Base: @propagate_inbounds
using Interpolations
using TimeSeries
using Dates
using Unitful
using Lazy
using DataStructures: SortedDict
using IntervalSets
using AxisArrays
using StaticArrays
using ComponentArrays
using Parameters
using FastClosures
using DifferentialEquations, DiffEqBase

# main version
include("core/core.jl")
include("layers/layers.jl")
include("processes/processes.jl")

# implicit version
# include("implicit/CryoGridImplicit.jl")

end # module
