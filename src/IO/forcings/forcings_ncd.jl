# NetCDF forcings formats
"""
    NCDFormat

Base type for NetCDF based forcing formats.
"""
abstract type NCDFormat end
"""
    NCDGeneric <: NCDFormat

Generic NetCDF forcing format.
"""
Base.@kwdef struct NCDGeneric <: NCDFormat
    timedim::String = "time"
end
"""
    NCDTopoPyScale

TopoPyScale NetCDF forcing format.
"""
struct NCDTopoPyScale <: NCDFormat end

"""
NetCDF forcing input format.
"""
Base.@kwdef struct ForcingFormatNCD{TFormat} <: ForcingFormat
    format::TFormat = NCDGeneric()
end
filesuffix(::ForcingFormatNCD) = "nc"

function detectformat(::Val{:nc}, filepath::String)
    NCD.Dataset(filepath) do ds
        source = get(ds.attrib, "source", "")
        if contains(source, "TopoPyScale")
            ForcingFormatNCD(NCDTopoPyScale())
        else
            ForcingFormatNCD()
        end
    end
end

function loadforcings(::ForcingFormatNCD{NCDTopoPyScale}, filepath::String)
    NCD.Dataset(filepath) do ds
        @assert haskey(ds, "time") "no time dimension found in NetCDF file"
        # global metadata
        title = get(ds.attrib, "title", "")
        source = get(ds.attrib, "source", "")
        creator_name = get(ds.attrib, "creator_name", "")
        date_created = get(ds.attrib, "date_created", "")
        # location data
        latitude = ds["latitude"][1]
        longitude = ds["longitude"][1]
        elevation = ds["elevation"][1]*uparse(ds["elevation"].attrib["units"])
        # time dimension
        timevar = ds["time"]
        calendar = timevar.attrib["calendar"]
        timestamps = collect(timevar)
        # select all variables with time as a dimension
        forcing_varnames = filter(name -> name != "time" && "time" ∈ NCD.dimnames(ds[name]), keys(ds))
        # map over variables names and construct forcing types from data
        forcings = map(forcing_varnames) do name
            var = ds[name]
            units = if haskey(var.attrib, "units")
                # parse units from string, assumes that all unit strings match Unitful.jl conventions;
                uparse(var.attrib["units"])
            else
                NoUnits
            end
            # add_offset = get(var.attrib, "add_offset", zero(eltype(var)))
            # scale_factor = get(var.attrib, "scale_factor", one(eltype(var)))
            # long_name = get(var.attrib, "long_name", name)
            standard_name = get(var.attrib, "standard_name", name)
            # note that this is loading the data entirely into memory;
            # may need to consider in the future lazy loading via DiskArray and/or the NetCDF package.
            # NetCDF.jl currently seems to have a bug where the data types are not inferred correctly.
            Symbol(name) => DimArray(var.*units, (Ti(timestamps),); name=Symbol(standard_name))
        end
        metadata = (; title, source, creator_name, date_created, latitude, longitude, elevation, calendar)
        # note the (;x...) syntax creates a named tuple from a collection of pairs, Symbol => value
        return DimStack((;forcings...); metadata)
    end
end
