"""
JSON forcing input format (from CryoGridLite) with specified version indicator.
"""
struct ForcingFormatJSON{Version} <: ForcingFormat end

filesuffix(::ForcingFormatJSON) = "json"

detectformat(::Val{:json}, filepath::String) = ForcingFormatJSON{1}()

forcingunits(::ForcingFormatJSON) = Dict(
    :Tair => u"°C",
    :pressure => u"Pa",
    :wind => u"m/s",
    :Lin => u"W/m^2",
    :Sin => u"W/m^2",
    :snowfall => u"mm/d",
    :rainfall => u"mm/d",
    :Ptot => u"mm/d",
)

function loadforcings(format::ForcingFormatJSON{1}, filename::String)
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
function loadforcings(format::ForcingFormatJSON{2}, filename::String)
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
