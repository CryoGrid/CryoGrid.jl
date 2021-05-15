module CryoGrid

using Base: @propagate_inbounds
using Lazy: @>>, @>, groupby
using DataStructures: SortedDict, OrderedDict
using Interpolations
using TimeSeries
using ArraysOfArrays
using StaticArrays
using ComponentArrays
using LinearAlgebra
using IfElse
using LoopVectorization
using RuntimeGeneratedFunctions
using Reexport
import ExprTools
import ForwardDiff
import ReverseDiff
@reexport using Dates
@reexport using Unitful
@reexport using IntervalSets
@reexport using DimensionalData
@reexport using Parameters
@reexport using DifferentialEquations
@reexport using DiffEqBase
@reexport using SimulationLogs

RuntimeGeneratedFunctions.init(CryoGrid)

# Temporary fix for bug in ExprTools (see issue #14)
# TODO: Remove when fixed in official package.
function ExprTools.argument_names(m::Method)
    slot_syms = ExprTools.slot_names(m)
    arg_names = slot_syms[2:m.nargs]  # nargs includes 1 for self ref
    return arg_names
end

include("io/io.jl")
include("core/core.jl")
include("layers/layers.jl")
include("processes/processes.jl")
include("models/Models.jl")

# Mandatory model method stubs
variables(::Layer) = ()
variables(::Layer, ::Process) = ()
initialcondition!(::Layer, state) = nothing
initialcondition!(::Layer, ::Process, state) = nothing
diagnosticstep!(l::Layer, p::Process, state) = error("no diagnostic step defined for $(typeof(l)) with $(typeof(p))")
prognosticstep!(l::Layer, p::Process, state) = error("no prognostic step defined for $(typeof(l)) with $(typeof(p))")
interact!(l1::Layer, p1::Process, l2::Layer, p2::Process, state1, state2) = nothing

export variables, parameters, initialcondition!, diagnosticstep!, prognosticstep!, interact!

end # module
