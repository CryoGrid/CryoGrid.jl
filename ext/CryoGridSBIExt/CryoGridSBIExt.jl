module CryoGridSBIExt

using ..SimulationBasedInference

using CryoGrid

# external packages
using Dates
using IntervalSets
using ModelParameters

# resolve ambiguous import warning
using ModelParameters: params

import Random

SimulationBasedInference.default_time_converter(::CryoGridProblem) = CryoGrid.convert_t

include("utils.jl")

export TemperatureProfileObservable, ActiveLayerThicknessObservable, LayerVarObservable
include("observables.jl")

export CryoGridPrior
include("priors.jl")

end
