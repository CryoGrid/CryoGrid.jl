"""
    loadforcings(filename::String, units...; spec::Type{T}=ForcingJSON{1})
    loadforcings(resource::Resource, units...; outdir=DEFAULT_FORCINGS_DIR)
    loadforcings(::Type{ForcingJSON{1}}, filename::String, units::Pair{Symbol,<:Unitful.Units}...)

Loads forcing data from the given file according to the format specified by `spec::ForcingInputFormat`. Default is ForcingJSON{1}.
Returns a NamedTuple of the form `(data=(...), timestamps=Array{DateTime,1})` where `data` is a NamedTuple matching
the structure of the JSON file. Units can (and should) be supplied as additional pair arguments, e.g:

`loadforcings("example.json", :Tair=>u"°C", :Ptot=>u"mm")`
"""
loadforcings(filename::String, units...; spec::Type{T}=_autodetect_forcing_format(filename)) where {T <: ForcingInputFormat} = loadforcings(T, filename, units...)
loadforcings(resource::Resource, units...; outdir=DEFAULT_FORCINGS_DIR) where {T <: ForcingInputFormat} = loadforcings(resource.format, fetch(resource, outdir), units...)
function loadforcings(::Type{ForcingJSON{1}}, filename::String, units::Pair{Symbol,<:Unitful.Units}...)
    dict = open(filename, "r") do file; JSON3.read(file) end
    # convert JSON3 dict for data field to Julia dict
    data = Dict(dict[:data]...)
    # get timestamps and then remove from dict
    ts = map(todatetime, map(DateNumber, data[:t_span]))
    delete!(data, :t_span)
    unitdict = Dict(units...)
    vals_with_units = (haskey(unitdict, name) ? replace(vals, nothing => missing)*unitdict[name] : vals for (name, vals) in data)
    # construct new named tuple
    (data = NamedTuple{Tuple(keys(data))}(tuple(vals_with_units...)), timestamps = Array(ts))
end
function loadforcings(::Type{ForcingJSON{2}}, filename::String, units::Pair{Symbol,<:Unitful.Units}...)
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
