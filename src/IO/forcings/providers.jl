"""
    ForcingProvider

Base type for forcing "providers", i.e. functional interfaces for
discrete forcing datasets.
"""
abstract type ForcingProvider <: InputProvider end

struct Interpolated1D{names,TVals} <: ForcingProvider
    data::DimStack
    interpolants::NamedTuple{names,TVals}
end

function Interpolated1D(
    forcingdata::DimStack,
    default_interp_mode=Linear();
    stripunits=true,
    interp_modes...)
    ts = convert_t.(dims(forcingdata, Ti))
    interpolants = map(forcingdata) do data
        @assert length(size(data)) == 1 "Interpolated1D only supports forcing data with a single dimension"
        data = stripunits ? ustrip.(data) : data
        interp = get(interp_modes, data.name, default_interp_mode)
        interpolate((ts,), data, Gridded(interp))
    end
    return Interpolated1D(forcingdata, interpolants)
end

(data::Interpolated1D)(input::Input, t::DateTime) = data(input, convert_t(t))
function (data::Interpolated1D)(::Input{name}, t::Number) where {name}
    itp = getproperty(data, name)
    return itp(t)
end

# special overrides for prroperty names and getproperty
Base.propertynames(::Interpolated1D{names}) = (:data, names...)
Base.getproperty(data::Interpolated1D, name::Symbol) = name == :data ? getfield(data, :data) : getproperty(inputmap(data), name)

function Base.show(io::IO, mimetype::MIME"text/plain", data::Interpolated1D)
    println("Interpolated1D of")
    show(io, mimetype, getfield(data, :data))
end

function inputmap(data::Interpolated1D)
    return map(_wrap_interpolant, getfield(data, :interpolants))
end

# wraps the interpolant to support both numeric and DateTime types.
function _wrap_interpolant(itp::Interpolations.AbstractInterpolation)
    interpolant(t::Number) = itp(t)
    interpolant(t::DateTime) = itp(convert_t(t))
    return interpolant
end
