module InputOutput

import CryoGrid

using CryoGrid.Numerics
using CryoGrid.Utils

using Base: @propagate_inbounds
using ComponentArrays
using ConstructionBase
using DataStructures: DefaultDict
using Dates
using DimensionalData
using Downloads
using Flatten
using Interpolations
using JSON3
using ModelParameters
using Tables
using Unitful

import DimensionalData

const INPUT_DIR = "input/"
const DEFAULT_FORCINGS_DIR = joinpath(INPUT_DIR, "forcings")
const DEFAULT_PARA_DIR = joinpath(INPUT_DIR, "para")

"""
    Resource{T}

Simple representation of a local or remote resource with a `name`, `type`, and `url`.
The `name` and `type` will be used to specify a file name and suffix respectively.
"""
struct Resource{T}
	name::String
	format::T
	url::String
end

"""
    fetch(resource::Resource, dir::String)

Downloads the remote `resource` from the location given by its URL to directory `dir`.
"""
function fetch(resource::Resource, dir::String)
	mkpath(dir)
	filepath = joinpath(dir, resource.name * "." * filesuffix(resource.format))
    if !isfile(filepath)
        Downloads.download(resource.url, filepath)
    else
        filepath
    end
end

include("ioutils.jl")
export CryoGridParams
include("params/params.jl")
export ParamsJSON, ParamsYAML
include("params/params_loaders.jl")
export Forcings, Forcing, ForcingJSON, ForcingNCD, ConstantForcing, InterpolatedForcing
export TemperatureForcing, WindForcing, HumidityForcing, EnergyFluxForcing, PressureForcing, VelocityForcing # aliases
include("forcings/forcings.jl")
export loadforcings
include("forcings/forcings_loaders.jl")
export CryoGridOutput
include("output.jl")

end
