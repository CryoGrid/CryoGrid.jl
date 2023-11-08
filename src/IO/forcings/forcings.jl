"""
      Forcing{unit,T}

Abstract type representing a generic external forcing term.
"""
abstract type Forcing{unit,T} end
Base.nameof(f::Forcing) = f.name
CryoGrid.parameterize(f::Forcing; ignored...) = f
@propagate_inbounds (forcing::Forcing)(x::Number) = error("$(typeof(forcing)) not implemented")
@propagate_inbounds (forcing::Forcing)(t::DateTime) = forcing(convert_t(t))

ConstructionBase.constructorof(::Type{F}) where {unit,F<:Forcing{unit}} = (args...) -> F.name.wrapper(unit, args...)

"""
Represents an externally specified format for forcing inputs. IO functions should dispatch on
specific types `T<:ForcingFormat` that they implement.
"""
abstract type ForcingFormat end

# Aliases for forcing types
const TemperatureForcing{T} = Forcing{u"°C",T} where {T}
const VelocityForcing{T} = Forcing{u"m/s",T} where {T}
const HeightForcing{T} = Forcing{u"m",T} where {T}
const HumidityForcing{T} = Forcing{u"kg/kg",T} where {T}
const PressureForcing{T} = Forcing{upreferred(u"Pa"),T} where {T}
const EnergyFluxForcing{T} = Forcing{upreferred(u"W/m^2"),T} where {T}

struct ConstantForcing{unit,T} <: Forcing{unit,T}
      value::T
      name::Symbol
      ConstantForcing(unit::Unitful.Units, value::T, name::Symbol) where {T} = new{unit,T}(value, name)
      ConstantForcing(qty::Unitful.AbstractQuantity, name::Symbol) = new{unit(qty),typeof(qty)}(qty, name)
      ConstantForcing(qty::Number, name::Symbol) = new{Unitful.NoUnits,typeof(qty)}(qty, name)
end
(f::ConstantForcing)(::Number) = f.value

"""
      InterpolatedForcing{unit,T,TI}

Forcing data provided by a discrete time series of data.
"""
struct InterpolatedForcing{unit,T,TI<:Interpolations.AbstractInterpolation{T}} <: Forcing{unit,T}
      interpolant::TI
      name::Symbol
      InterpolatedForcing(unit::Unitful.Units, interpolant::TI, name::Symbol) where {TI} = new{unit,eltype(interpolant),TI}(interpolant, name)
end
function InterpolatedForcing(timestamps::AbstractArray{DateTime,1}, values::A, name::Symbol; interpolation_mode=Interpolations.Linear()) where {T,A<:AbstractArray{T,1}}
      ts = convert_t.(timestamps)
      values_converted = Utils.normalize_units.(values)
      interp = Interpolations.interpolate((ts,), ustrip.(values_converted), Interpolations.Gridded(interpolation_mode))
      return InterpolatedForcing(unit(eltype(values_converted)), interp, name)
end
Flatten.flattenable(::Type{<:InterpolatedForcing}, ::Type) = false
Base.getindex(f::InterpolatedForcing, i::Integer) = f.interpolant.coefs[i]
Base.getindex(f::InterpolatedForcing, t) = f(t)
Base.length(f::InterpolatedForcing) = length(f.interpolant)
Base.show(io::IO, forcing::InterpolatedForcing{u}) where u = print(io, "InterpolatedForcing $(forcing.name) [$u] of length $(length(forcing)) with time span $(convert_tspan(extrema(forcing.interpolant.knots[1])))")
# Base.show(io::IO, ::Type{<:InterpolatedForcing}) = print(io, "InterpolatedForcing") 
# Base.show(io::IO, ::Type{<:InterpolatedForcing{unit}}) where {unit} = print(io, "InterpolatedForcing{$unit}") 
# Base.show(io::IO, ::Type{<:InterpolatedForcing{unit,T}}) where {unit,T} = print(io, "InterpolatedForcing{$unit, $T, ...}") # type piracy to truncate type name in stacktraces
"""
Get interpolated forcing value at t seconds from t0.
"""
@propagate_inbounds (forcing::InterpolatedForcing)(t::Number) = forcing.interpolant(t)

"""
    time_derivative_forcing(
        f::InterpolatedForcing{unit},
        new_name::Symbol;
        interp=Numerics.Linear()
    ) where {unit}

Computes the finite difference time derivative of the given `InterpolatedForcing`
time series and returns a new forcing with units `[unit].s⁻¹`
"""
function time_derivative_forcing(
    f::InterpolatedForcing{unit},
    new_name::Symbol;
    interp=Numerics.Linear()
) where {unit}
    ts = f.interpolant.knots
    ∂f∂t = map(t -> Numerics.gradient(f, t)*unit/1.0u"s", ts)
    ts_mid = (ts[1:end-1] .+ ts[2:end])./2
    return InterpolatedForcing(convert_t.(ts_mid), ∂f∂t, new_name; interp)
end

"""
      Forcings{names,TF,TMeta}

Generic container for forcing types with optional metadata.
"""
struct Forcings{names,TF,TMeta}
      data::NamedTuple{names,TF}
      metadata::TMeta
      Forcings(data::NamedTuple{names,TF}, metadata::NamedTuple=(;)) where {names,TF} = new{names,TF,typeof(metadata)}(data, metadata)
end
Flatten.flattenable(::Type{<:Forcings}, ::Val{:metadata}) = false
CryoGrid.parameterize(f::Forcings; ignored...) = f
Base.propertynames(::Forcings{names}) where {names} = (:metadata, names...)
Base.getproperty(fs::Forcings{names}, sym::Symbol) where {names} = sym ∈ names ? getproperty(getfield(fs, :data), sym) : getfield(fs, sym)
Base.merge(fs::Forcings...) = Forcings(merge(map(f -> getfield(f, :data), fs)...), merge(map(f -> getfield(f, :metadata), fs)...))

forcingunits(::ForcingFormat) = Dict()

detectformat(::Val{x}, filepath) where x = error("unrecognized forcing file suffix $x")
function detectformat(filepath::String)
    filename = basename(filepath)
    suffix = Symbol(lowercase(split(filename, ".")[end]))
    return detectformat(Val{suffix}(), filepath)
end

"""
    loadforcings(filename::String)::Forcings
    loadforcings(resource::Resource; outdir=DEFAULT_FORCINGS_DIR)::Forcings
    loadforcings([format::ForcingFormat], filename::String; outdir=DEFAULT_FORCINGS_DIR)::Forcings

Loads forcing data from the given file according to the format specified by `format`. By default, the forcing format
is automatically detected via `detectformat`. Returns a `Forcings` struct containing all forcing data
and metadata 
"""
loadforcings(filename::String) = loadforcings(detectformat(filename), filename)
loadforcings(resource::Resource; outdir=DEFAULT_FORCINGS_DIR) = loadforcings(resource.format, fetch(resource, outdir))
loadforcings(f::ForcingFormat, filename::String) = error("loadforcings not implemented for $f")

_normalize_numeric(x::Number) = convert(Float64, x)
_normalize_numeric(::Union{Missing,Nothing}) = missing

include("forcings_json.jl")
include("forcings_ncd.jl")
