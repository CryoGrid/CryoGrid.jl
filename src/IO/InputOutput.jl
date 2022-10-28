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
using Lazy: @>>, groupby
using ModelParameters
using Tables
using TimeSeries
using Unitful

import DimensionalData

export ForcingJSON, ParamJSON, ParamYAML

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

"""
Represents an externally specified format for forcing inputs. IO functions should dispatch on
specific types `T<:ForcingInputFormat` that they implement.
"""
abstract type ForcingInputFormat end
"""
JSON forcing input format (from CryoGridLite). Version parameter allows for forward/backward-compatibility with future format changes.
"""
struct ForcingJSON{Version} <: ForcingInputFormat end

abstract type ParameterInputFormat end
"""
JSON parameter input format (from CryoGridLite).
"""
struct ParamJSON{Version} <: ParameterInputFormat end
"""
YAML parameter input format matching that of the CryoGrid community model.
"""
struct ParamYAML{Version} <: ParameterInputFormat end

filesuffix(::Type{<:ForcingJSON}) = "json"
filesuffix(::Type{<:ParamJSON}) = "json"
filesuffix(::Type{<:ParamYAML}) = "yml"

function _autodetect_forcing_format(filepath::String)
    filename = basename(filepath)
    if endswith(filename, ".json") || endswith(filename, ".JSON")
        return ForcingJSON{1}
    else
        error("unsupported forcing file format: $filename")
    end
end

include("ioutils.jl")
export CryoGridParams
include("params.jl")
export Forcing, TimeSeriesForcing
include("forcings.jl")
export CryoGridOutput
include("output.jl")
export loadforcings
include("forcings_loaders.jl")
include("params_loaders.jl")

end
