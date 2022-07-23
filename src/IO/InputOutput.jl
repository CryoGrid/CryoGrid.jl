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

import DimensionalData: stack

export loadforcings

include("ioutils.jl")
export CryoGridParams, parameterize
include("params.jl")
export Forcing, TimeSeriesForcing
include("forcings.jl")
export CryoGridOutput
include("output.jl")

const INPUT_DIR = "input/"
const DEFAULT_FORCINGS_DIR = joinpath(INPUT_DIR, "forcings")
const DEFAULT_PARA_DIR = joinpath(INPUT_DIR, "para")

struct Resource
	name::String
	type::String
	url::String
end

function fetch(resource::Resource, dir::String)
	mkpath(dir)
	filepath = joinpath(dir, resource.name * "." * resource.type)
    if !isfile(filepath)
        Downloads.download(resource.url, filepath)
    else
        filepath
    end
end

"""
Represents an externally specified format for inputs (i.e. forcings and parameters). IO functions should dispatch on
specific types `T<:InputSpec` that they implement.
"""
abstract type InputSpec end
"""
JSON input format. Version parameter allows for forward/backward-compatibility with future format changes.
"""
struct JsonSpec{Version} <: InputSpec end

export InputSpec, JsonSpec

"""
    loadforcings(filename::String, units...; spec::Type{T}=JsonSpec{1})
    loadforcings(resource::Resource, units...; spec::Type{T}=JsonSpec{1}, outdir=DEFAULT_FORCINGS_DIR)
    loadforcings(::Type{JsonSpec{1}}, filename::String, units::Pair{Symbol,<:Unitful.Units}...)

Loads forcing data from the given file according to the format specified by `spec::InputSpec`. Default is JsonSpec{1}.
Returns a NamedTuple of the form `(data=(...), timestamps=Array{DateTime,1})` where `data` is a NamedTuple matching
the structure of the JSON file. Units can (and should) be supplied as additional pair arguments, e.g:

`loadforcings("example.json", :Tair=>u"°C", :Ptot=>u"mm")`
"""
loadforcings(filename::String, units...; spec::Type{T}=JsonSpec{1}) where {T <: InputSpec} = loadforcings(T, filename, units...)
loadforcings(resource::Resource, units...; spec::Type{T}=JsonSpec{1}, outdir=DEFAULT_FORCINGS_DIR) where {T <: InputSpec} = loadforcings(T, fetch(resource, outdir), units...)
function loadforcings(::Type{JsonSpec{1}}, filename::String, units::Pair{Symbol,<:Unitful.Units}...)
    dict = open(filename, "r") do file; JSON3.read(file) end
    # convert JSON3 dict for data field to Julia dict
    data = Dict(dict[:data]...)
    # get timestamps and then remove from dict
    ts = @>> data[:t_span] map(DateNumber) map(todatetime)
    delete!(data, :t_span)
    unitdict = Dict(units...)
    vals_with_units = (haskey(unitdict, name) ? vals*unitdict[name] : vals for (name, vals) in data)
    # construct new named tuple
    (data = NamedTuple{Tuple(keys(data))}(tuple(vals_with_units...)), timestamps = Array(ts))
end
function loadforcings(::Type{JsonSpec{2}}, filename::String, units::Pair{Symbol,<:Unitful.Units}...)
    dict = open(filename, "r") do file; JSON3.read(file) end
    # convert JSON3 dict for data field to Julia dict
    data = Dict(dict[:data]...)
    # get timestamps
    ts = map(DateTime, dict[:timestamps])
    unitdict = Dict(units...)
    vals_with_units = (haskey(unitdict, name) ? vals*unitdict[name] : vals for (name, vals) in data)
    # construct new named tuple
    (data = NamedTuple{Tuple(keys(data))}(tuple(vals_with_units...)), timestamps = Array(ts))
end

end
