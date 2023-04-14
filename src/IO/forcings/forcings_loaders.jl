"""
Represents an externally specified format for forcing inputs. IO functions should dispatch on
specific types `T<:ForcingDataFormat` that they implement.
"""
abstract type ForcingDataFormat end
"""
JSON forcing input format (from CryoGridLite) with specified version indicator.
"""
struct ForcingJSON{Version} <: ForcingDataFormat end
"""
NetCDF forcing input format.
"""
struct ForcingNCD{Version} <: ForcingDataFormat end

filesuffix(::ForcingJSON) = "json"

forcingunits(::ForcingDataFormat) = Dict()
forcingunits(::ForcingJSON) = Dict(
    :Tair => u"°C",
    :pressure => u"Pa",
    :wind => u"m/s",
    :Lin => u"W/m^2",
    :Sin => u"W/m^2",
    :snowfall => u"mm/d",
    :rainfall => u"mm/d",
    :Ptot => u"mm",
)

function autodetect_forcing_format(filepath::String)
    filename = basename(filepath)
    if endswith(filename, ".json") || endswith(filename, ".JSON")
        return ForcingJSON{1}()
    else
        error("unsupported forcing file format: $filename")
    end
end

_normalize_numeric(x::Number) = convert(Float64, x)
_normalize_numeric(::Union{Missing,Nothing}) = missing

"""
    loadforcings(filename::String)::Forcings
    loadforcings(resource::Resource; outdir=DEFAULT_FORCINGS_DIR)::Forcings
    loadforcings([format::ForcingDataFormat], filename::String; outdir=DEFAULT_FORCINGS_DIR)::Forcings

Loads forcing data from the given file according to the format specified by `format`. By default, the forcing format
is automatically detected via `autodetect_forcing_format`. Returns a `Forcings` struct containing all forcing data
and metadata 
"""
loadforcings(filename::String) = loadforcings(autodetect_forcing_format(filename), filename)
loadforcings(resource::Resource; outdir=DEFAULT_FORCINGS_DIR) = loadforcings(resource.format, fetch(resource, outdir))
function loadforcings(format::ForcingJSON{1}, filename::String, units::Pair{Symbol,<:Unitful.Units}...)
    dict = open(filename, "r") do file; JSON3.read(file) end
    # convert JSON3 dict for data field to Julia dict
    data = Dict(dict[:data]...)
    names = keys(data)
    # get timestamps and then remove from dict
    ts = map(todatetime, map(DateNumber, data[:t_span]))
    delete!(data, :t_span)
    unitdict = forcingunits(format)
    # due to some irregularities in the forcing format, we need to replace 'nothing' values with 'missing'
    vals_with_units = map(names, values(data)) do name, vals
        vals = map(_normalize_numeric, vals)
        if haskey(unitdict, name)
            vals = vals*unitdict[name]
        else
            vals
        end
    end
    forcings = map(names, vals_with_units) do name, values
        name => InterpolatedForcing(Array(ts), values, name)
    end
    return Forcings((; forcings...))
end
function loadforcings(format::ForcingJSON{2}, filename::String, units::Pair{Symbol,<:Unitful.Units}...)
    dict = open(filename, "r") do file; JSON3.read(file) end
    # convert JSON3 dict for data field to Julia dict
    data = Dict(dict[:data]...)
    names = keys(data)
    # get timestamps
    ts = map(DateTime, dict[:timestamps])
    unitdict = forcingunits(format)
    vals_with_units = map(names, values(data)) do name, vals
        vals = map(_normalize_numeric, vals)
        if haskey(unitdict, name)
            vals = vals*unitdict[name]
        else
            vals
        end
    end
    forcings = map(names, vals_with_units) do name, values
        name => InterpolatedForcing(Array(ts), values, name)
    end
    return Forcings((; forcings...))
end
