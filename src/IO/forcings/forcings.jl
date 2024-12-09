"""
Represents an externally specified format for forcing inputs. IO functions should dispatch on
specific types `T<:ForcingFormat` that they implement.
"""
abstract type ForcingFormat end

forcingunits(::ForcingFormat) = Dict()

detectformat(::Val{x}, filepath) where x = error("unrecognized forcing file suffix $x")
function detectformat(filepath::String)
    filename = basename(filepath)
    suffix = Symbol(lowercase(split(filename, ".")[end]))
    return detectformat(Val{suffix}(), filepath)
end

"""
    loadforcings(filename::String)
    loadforcings(resource::Resource; outdir=DEFAULT_FORCINGS_DIR)
    loadforcings([format::ForcingFormat], filename::String; outdir=DEFAULT_FORCINGS_DIR)

Loads forcing data from the given file according to the format specified by `format`. By default, the forcing format
is automatically detected via `detectformat`. Returns a `DimStack` containing all forcing data
and metadata 
"""
loadforcings(filename::String) = loadforcings(detectformat(filename), filename)
loadforcings(resource::Resource; outdir=DEFAULT_FORCINGS_DIR) = loadforcings(resource.format, fetch(resource, outdir))
loadforcings(f::ForcingFormat, filename::String) = error("loadforcings not implemented for $f")

_normalize_numeric(x::Number) = convert(Float64, x)
_normalize_numeric(::Union{Missing,Nothing}) = missing

include("providers.jl")
include("forcings_json.jl")
include("forcings_ncd.jl")
