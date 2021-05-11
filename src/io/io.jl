using JSON3
using DataStructures: DefaultDict

include("ioutils.jl")

const INPUT_DIR = "input/"
const FORCINGS_DIR = joinpath(INPUT_DIR, "forcings")
const PARA_DIR = joinpath(INPUT_DIR, "para")

struct Resource
	name::String
	type::String
	url::String
end

function fetch(resource::Resource, dir::String)
	mkpath(dir)
	filepath = joinpath(dir, resource.name*"."*resource.type)
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
Loads forcing data from the given file according to the format specified by `spec::InputSpec`. Default is JsonSpec{1}.
Returns a NamedTuple of the form `(data=(...), timestamps=Array{DateTime,1})` where `data` is a NamedTuple matching
the structure of the JSON file. Units can (and should) be supplied as additional pair arguments, e.g:

`loadforcings("example.json", :Tair=>u"°C", :Ptot=>u"mm")`
"""
loadforcings(filename::String, units...; spec::Type{T}=JsonSpec{1}) where {T<:InputSpec} = loadforcings(T, filename, units...)
loadforcings(resource::Resource, units...; spec::Type{T}=JsonSpec{1}) where {T<:InputSpec} = loadforcings(T, fetch(resource, FORCINGS_DIR), units...)
function loadforcings(::Type{JsonSpec{1}}, filename::String, units::Pair{Symbol,<:Unitful.Units}...)
        dict = open(filename,"r") do file; JSON3.read(file) end
        # convert JSON3 dict for data field to Julia dict
        data = Dict(dict[:data]...)
        # get timestamps and then remove from dict
        ts = @>> data[:t_span] map(DateNumber) map(todatetime)
        delete!(data, :t_span)
        unitdict = Dict(units...)
        # Case 1: Map number values to units, if specified
        tounit(val::Number,name) = val*(haskey(unitdict,name) ? unitdict[name] : unit(0))
        # Case 2: Map missing values to missing
        tounit(val::Nothing,name) = missing
        unitvals = (tounit.(vals,name) for (name,vals) in data)
        # construct new named tuple
        (data=NamedTuple{Tuple(keys(data))}(tuple(unitvals...)), timestamps=Array(ts))
end

export loadforcings
